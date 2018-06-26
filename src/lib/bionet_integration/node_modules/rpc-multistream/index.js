var through = require('through2');
var isStream = require('isstream');
var pump = require('pump');
var Multiplex = require('multiplex');
var uuid = require('uuid').v4;
var xtend = require('xtend');
var async = require('async');
var jsonStream = require('duplex-json-stream');

/*
  Below is some documentation for the message format used for RPC calls. 
  An even better way to learn this is to set debug to true which will
  print every message to stdout. Like so:

    var rpc = require('rpc-multistream');
    var myRpc = rpc({ ... }, {debug: true});

  rpcStream message format:
  
  ['functionNameOrCallbackId', [functionArgs], {indexOfCallbackArg1: 'callback1Id', indexOfCallbackArg2: 'callback2Id'}, {indexOfStreamArg1: stream1Identifier, indexOfStreamArg2: stream2Identifier}, returnStreamIdenfifier]
 
  Everything after functionNameOrCallbackId is optional but position in the array is important. So you can do [val1, val2] and leave out the rest but not [val1, val2, val4] without doing [val1, val2, undefined, val4].

  returnStreamIdentifier can also be an array if you plan to return multiple streams using synchronous calling.

  A streamIdentifier is of the following format:

  [streamId, 'type', objectMode, 'encoding']

  where 'type' is 'd', 'r', or 'w'
  and can be skipped, in which case type defaults to 'd' duplex:

  ['streamId', objectMode, 'encoding']

  or just

  ['streamId']

  and for all undefined stream options, 
  the values from the rpc-multiplex opts are used.

  Example: 

  A method is called like so:
  
  var myStream = fs.createReadStream('/dev/urandom');
  remoteMethod.foo('hello', myStream, function(err) {
  // ...
  });
  
  Both the callback function and myStream will be saved locally 
  and assigned a unique id. 
  Let's say myStream gets id 101 and the callback id 202.
  The message sent will look like this:
  
  ['foo', ['hello'], {1: [101}, {2: 202}]
  
  If foo returned a stream, then that stream would also be given an id,
  say it got id 303, the message would look like this:

  ['foo', ['hello'], {1: 101}, {2: 202}, 303]
  
  or foo could return multiple streams, in which case we'd see something like:
  
  ['foo', ['hello'], {1: 101}, {2: 202}, [303, 404, 505]]
  
  Once the callback got called, we'd get back a message like:
  
  ['202', [null, 'remote says hi']]
  
  Though there is no reason why a callback couldn't also include streams
  both as arguments and as return values.

*/

var moduleName = 'rpc-multistream';

function streamType(stream) {
  if(!isStream(stream)) return null;
  if(isStream.isDuplex(stream)) return 'd';
  if(isStream.isWritable(stream)) return 'w';
  if(isStream.isReadable(stream)) return 'r';
  return null;
}

function flattenError(err) {
  if(!(err instanceof Error)) return err;
  var err2 = {
    message: err.message
  };
  err2._isErrorObject = true;

  for(var k in err) {
    err2[k] = err[k] 
  }
  return err2
};

function expandError(err) {
  if (!err || !err._isErrorObject) return err
  var err2 = new Error(err.message);
  for(var k in err) {
    if(k === '_isErrorObject') continue;
    err2[k] = err[k]
  }
  return err2;
};


function rpcMultiStream(methods, opts) {
  opts = xtend({
    init: true, // automatically send rpc methods manifest on instantiation
    // TODO implement detectEncoding
    detectEncoding: true, // detect encoding and objectMode for streams
    // TODO implement detectStreamType
    detectStreamType: true, // detect if streams are readable/writeable/duplex
    encoding: 'utf8', // default encoding for streams
    objectMode: false, // default objectMode for streams
    explicit: true, // include encoding/objectMode even if they match defaults
    debug: false,
    heartbeat: 0, // send heartbeat every n ms. disable if falsy
    maxMissedBeats: 3.5, // die after missing this many heartbeats
    flattenError: flattenError,
    expandError: expandError,
    onError: function(err, functionName, multiplex, metaStream) {
      // emit error both on local and remote sides
      multiplex.emit('error', err, false, functionName);
      if(typeof opts.flattenError === 'function') {
        err = opts.flattenError(err);
      }
      metaStream.write(['error', {err: err, functionName: functionName}]);
    }
  }, (opts || {}));

  var diceRoll = uuid();
  var isEven;

  var methods = methods;

  var multiplex = Multiplex(opts);

  var callbacks = []; // saved callbacks
  var cbCount;

  var streams = []; // saved streams
  var streamCount; // set when manifest received

  var metaStream = makeStream('m', {objectMode: true});
  var rpcStream = makeStream('r', {objectMode: true});
  var heartbeatStream = makeStream('h');

  var playingDead = false;
  var lastHeartbeat;
  heartbeatStream.on('data', function(data) {
    if(playingDead) return;
    var found42 = false;
    var found43 = !opts.heartbeat;
    var i;
    for(i=0; i < data.length; i++) {
      if(found42 && found43) break;
      if(!found42 && data[i] === 42) {
        heartbeatStream.write(new Buffer([43]))
        found42 = true;
        continue;
      }
      if(!found43 && data[i] === 43) {
        lastHeartbeat = new Date().getTime();
        multiplex.emit('heartbeat', lastHeartbeat);
        found43 = true;
      }
    }
  });

  metaStream.on('data', handleMeta);

  metaStream.on('error', function(err) {
    multiplex.emit('error', err);
  });

  rpcStream.on('data', function(data) {
    if(!(data instanceof Array) || data.length < 1) {
      return;
    }
    handleRPC(data);
  });

  rpcStream.on('error', function(err) {
    multiplex.emit('error', err);
  });

  function heartbeat() {
    if(!opts.heartbeat) return;
    if(lastHeartbeat && (new Date()).getTime() - lastHeartbeat > Math.ceil(opts.heartbeat * opts.maxMissedBeats)) {
      multiplex.emit('death');
      return;
    }
    heartbeatStream.write(new Buffer([42]));
    setTimeout(heartbeat, opts.heartbeat);
  }

  if(opts.heartbeat) {
    heartbeat();
  }


  if(opts.init) {
    init();
  }

  // --- functions below

  function inc(count) {
    count += 2;
    if(count >= 9007199254740991) {
      count = (isEven) ? 2 : 1; // this will almost certainly never happen
    }
    return count;
  }

  function restoreCallbackArgs(args, cbArgs) {
    var i;
    for(i in cbArgs) {
      args[i] = createRemoteCall(cbArgs[i]);
    }
  }

  function streamFromDesc(desc, functionName) {
    desc = desc.slice(0); // clone
    var id;
    var sopts = {
      encoding: opts.encoding,
      objectMode: opts.objectMode
    };
    var type = 'd';
    id = parseInt(desc.shift());
    if(typeof desc[0] === 'string') {
      type = desc.shift();
    }
    if(typeof desc[0] === 'boolean') {
      sopts.objectMode = desc[0];
    }
    if(typeof desc[1] !== 'undefined') {
      sopts.encoding = desc[1];
    }
    if(!id) throw new Error("Non-numeric stream id");
    var rs = makeStream(id, type, sopts);

    rs.on('error', function(err) {
      if(err.message.match(/^channel destroyed/i)) return;

      if(typeof opts.flattenError === 'function') {
        err = opts.flattenError(err);
      }
      metaStream.write(['streamError', {streamId: id, err: err}]);
    });

    return rs;
  }

  function restoreStreamArgs(args, streamDescs) {
    var i;
    for(i in streamDescs) {
      args[i] = streamFromDesc(streamDescs[i])
    }
  }

  function restoreRetStreams(streamDescs) {
    var ret = [];
    var i, desc;
    for(i=0; i < streamDescs.length; i++) {
      desc = streamDescs[i];
      ret.push(streamFromDesc(desc));
      streamCount = inc(streamCount);
    }
    return ret;
  }

  // errors that have no means of being reported
  // to the remote end will be sent here
  function uncaughtError(err, functionName) {
    if(typeof opts.onError === 'function') {
      return opts.onError(err, functionName, multiplex, metaStream);
    }
  }

  function debug() {
    if(!opts.debug) return;
    var args = [].slice.call(arguments);
    args = ['['+moduleName+' debug]'].concat(args);
    console.log.apply(console, args);
  }

  function handleRPC(data) {
    debug("receiving on rpcStream:", data);
    var fn;
    var cbi = parseInt(data[0]);
    if(cbi) {
      fn = callbacks[cbi];
    } else {
      var name = data[0];
      fn = methods[name];
    }
    if(!fn) return console.error("Invalid RPC:", data);

    var args = data[1] || [];
    var cbArgs = data[2];
    var streamArgs = data[3];
    var retStreamDescs = data[4];
    var retStreams;

    function handleError(err) {
      // if last argument is a function then assume it's the main callback
      // and call it with the error as only argument
      if(args.length && typeof args[args.length-1] === 'function') {
        if(typeof opts.flattenError === 'function') {
          err = opts.flattenError(err);
        }
        args[args.length-1].call(methods, err);

        // if we have at least one stream then we use that to report the error
      } else if(retStreams && retStreams.length) { 
        if(typeof opts.flattenError === 'function') {
          err = opts.flattenError(err);
        }
        var reported = false;
        var streamId;
        if(!retStreams) return; // TODO error here?

        for(i=0; i < retStreams.length; i++) {
          streamId = retStreamDescs[i][0];
          metaStream.write(['streamError', {streamId: streamId, functionName: name, err: err}]);
          // destroy streams that didn't get connected
          if(!ret || !ret[i]) {
            retStreams[i].destroy();
          }
        }
      } else {
        uncaughtError(err, name);
      }
    }

    try {
      // expand errors only for first argument to callback
      if(cbi) { 
        if(args.length && (typeof opts.expandError === 'function')) {
          args[0] = opts.expandError(args[0]);
        }
      }
      
      if(cbArgs) restoreCallbackArgs(args, cbArgs);
      if(streamArgs) restoreStreamArgs(args, streamArgs);

      if(retStreamDescs && retStreamDescs.length) {
        if(!(retStreamDescs[0] instanceof Array)) {
          retStreamDescs = [retStreamDescs];
        }
        retStreams = restoreRetStreams(retStreamDescs);
      }

      var ret = fn.apply(methods, args);
      if(!(ret instanceof Array)) ret = [ret];

      if(retStreams) {
        var type, rs;
        var i = 0;
        async.eachSeries(retStreams, function(retStream, cb) {
          try {

            rs = ret[i];
            if(!rs) return cb();

            rs.on('error', function(err) {
              retStream.emit('error', err);
            });

            type = retStreamDescs[i][1];
            if(type === 'r') {
              // TODO switch to pump
              rs.pipe(retStream);
              // pump(rs, retStream);
            } else if(type === 'w') {
              retStream.pipe(rs);

              // pump(retStream, rs);
            } else { // duplex
              rs.pipe(retStream).pipe(rs);
              // pump(rs, retStream, rs);
            }
            i++;
            cb();

          } catch(err) {
            handleError(err);
          }
        }, function(err) {
          if(err) return handleError(err);
        });
      }

    } catch(err) {
      handleError(err);
    }
  }

  /*
    same as multiplex.createSharedStream
    but supports the additional opts:
    * encoding
    * objectMode
    */
  function makeStream(id, type, opts) {
    if(typeof type === 'object') {
      opts = type;
      type = 'd';
    }
    opts = opts || {};

    // TODO close unused end of non-duplex streams
    // remember that type can be 'd' or 'duplex' etc.
    var stream = multiplex.createSharedStream(id, opts);

    if(opts.encoding) {
      stream.setEncoding(opts.encoding);
      stream.setDefaultEncoding(opts.encoding);
    }
    if(opts.objectMode) {
      var jstream = jsonStream(stream);

      // forward error
      jstream.on('error', function(err) {
        stream.emit('error', err);
      }); 
      return jstream;
    }

    return stream;
  };


  function genManifest(methods) {
    // each side rolls a 2^128 sided dice
    // whoever gets the higher number uses even indexes
    // the other uses uneven indexes
    var manifest = {
      '.diceRoll': diceRoll
    };
    var name;
    for(name in methods) {
      if(methods[name]._rpcOpts) {
        methods[name]._rpcOpts = xtend({
          encoding: opts.encoding,
          objectMode: opts.objectMode,
          type: 'd'
        }, methods[name]._rpcOpts);
        manifest[name] = methods[name]._rpcOpts;
      } else {
        manifest[name] = 'async';
      }
    }
    return manifest;
  }

  function streamOpts(stream) {
    var sopts = {};
    // TODO care about different opts for readable and writable?
    var state = stream._readableState || stream._writableState;
    if(state.objectMode) {
      sopts.objectMode = true;
      sopts.encoding = 'utf8';
    } else {
      sopts.encoding = ('encoding' in state) ? state.encoding : opts.encoding;
      sopts.objectMode = false
    }
    return sopts;
  }

  function streamIdentifier(streamId, type, sopts) {
    var id = [streamId];
    
    if(type !== 'd') {
      id.push(type);
    }

    if((sopts.encoding === opts.encoding) && !opts.explicit) {
      if(sopts.objectMode != opts.objectMode) {
        id.push(sopts.objectMode);
      }
    } else {
      id.push(sopts.objectMode);
      id.push(sopts.encoding);
    }

    return id;
  }
  
  function registerStreams(args) {
    var mapping = {};
    var i, s, type, opts, stream;
    for(i=0; i < args.length; i++) {
      s = args[i];
      type = streamType(s);
      if(!type) continue; // not a stream

      opts = streamOpts(s);
      stream = makeStream(streamCount, type, opts);
      streams[streamCount] = stream;
      mapping[i] = streamIdentifier(streamCount, type, opts);
      args[i] = null;
      streamCount = inc(streamCount);;
      if(type === 'r') {
        // TODO what to do on error here?
        pump(s, stream);
      } else if(type === 'w') {
        pump(stream, s);
      } else {
        pump(s, stream, s);
      }
    }
    if(!Object.keys(mapping).length) return null;
    return mapping;
  }

  function registerCallbacks(args) {
    var mapping = {};
    var i, cb;
    for(i=0; i < args.length; i++) {
      cb = args[i];
      if(typeof cb !== 'function') continue;

      callbacks[cbCount] = cb
      mapping[i] = cbCount;
      args[i] = null;
      cbCount = inc(cbCount);
    }
    if(!Object.keys(mapping).length) return null;
    return mapping;
  }

  // register streams used as return values
  function registerReturnStreams(streamOpts) {
    var ids = [];
    var retStreams = [];
    
    var i, type, sopts, stream;
    for(i=0; i < streamOpts.length; i++) {
      sopts = {
        objectMode: (typeof streamOpts[i].objectMode === 'boolean') ? streamOpts[i].objectMode : opts.objectMode,
        encoding: ('encoding' in streamOpts[i]) ? streamOpts[i].encoding : opts.encoding
      };
      type = streamOpts[i].type || 'd';

      ids.push(streamIdentifier(streamCount, type, sopts));

      stream = makeStream(streamCount, type, sopts);

      retStreams.push(stream);
      streams[streamCount] = stream;

      streamCount = inc(streamCount);
    }
    return {
      ids: (ids.length > 1) ? ids : ids[0],
      streams: (retStreams.length > 1) ? retStreams : retStreams[0]
    };
  }
  function createRemoteCall(name, retStreamOpts) {
    return function() {
      var args = [].slice.call(arguments);
      var cbMapping = registerCallbacks(args);
      var streamMapping = registerStreams(args);
      var returnMapping = null;
      
      var msg = [name];
      if(args && args.length) {
        // flatten error only for first argument of callbacks 
        if(parseInt(name) && typeof opts.flattenError === 'function') {
          args[0] = opts.flattenError(args[0]);
        }
      }
      msg.push(args);
      if(cbMapping) msg.push(cbMapping);
      if(streamMapping) {
        if(!cbMapping) msg.push({});
        msg.push(streamMapping);
      }
      if(retStreamOpts && retStreamOpts !== 'async') {
        if(!(retStreamOpts instanceof Array)) {
          retStreamOpts = [retStreamOpts];
        }
        if(!cbMapping) msg.push({});
        if(!streamMapping) msg.push({});

        var ret = registerReturnStreams(retStreamOpts);
        msg.push(ret.ids);
        debug("sending on rpcStream:", typeof msg, msg);
        rpcStream.write(msg);

        return ret.streams;
      }
      debug("sending on rpcStream:", typeof msg, msg);
      rpcStream.write(msg);
    }
  }

  function initIndexes(even) {
    isEven = even;
    if(even) {
      cbCount = 2;
      streamCount = 2;
    } else {
      cbCount = 1;
      streamCount = 1;
    }
  }

  function gotManifest(manifest) {
    if(!manifest['.diceRoll']) uncaughtError("Manifest is missing diceRoll");
    if(manifest['.diceRoll'] > diceRoll) {
      initIndexes(true);            
    } else if(manifest['.diceRoll'] < diceRoll) {
      initIndexes(false);
    } else {
      multiplex.destroy("Both endpoints generated the same roll on a 2^128 sided dice. Congratulations?");
      return;
    }

    var name;
    var methods = {};
    for(name in manifest) {
      if(name[0] === '.') continue;
      methods[name] = createRemoteCall(name, manifest[name]);
    }
    multiplex.emit('methods', methods);
  }
  
  function sendMeta(msgType, msg) {
    var out = [msgType, msg];
    debug("sending on metaStream:", out);
    metaStream.write(out);
  }

  function handleMeta(data) {
    debug("receiving on metaStream:", data);
    if(!(data instanceof Array)|| data.length < 2) return;
    var msgType = data[0];
    var msg = data[1];
    if(!msgType || !msg) return;

    switch(msgType) {
      
    case 'manifest':
      gotManifest(msg);
      break;
    case 'streamError': // error to be emitted on a stream
      var s = streams[msg.streamId];
      if(s) {
        if(typeof opts.expandError === 'function') {
          msg.err = opts.expandError(msg.err);
        }
        s.emit('error', msg.err, true, msg.functionName);
        break;
      }
    case 'error': // error that had no callback or stream for getting reported
      if(typeof opts.expandError === 'function') {
        msg.err = opts.expandError(msg.err);
      }
      multiplex.emit('error', msg.err, true, msg.functionName);
      break;
    default:
      multiplex.emit('error', "Unknown meta message type: " + msgType);
    }
  }

  function init(cb) {
    var manifest = genManifest(methods);
    if(manifest) {
      sendMeta('manifest', manifest);
    }
  };

  // restart heartbeat
  multiplex.revive = function(_opts) {
    _opts = _opts || {};
    if(_opts.heartbeat) opts.heartbeat = _opts.heartbeat;
    if(_opts.maxMissedBeats) opts.maxMissedBeats = _opts.maxMissedBeats;
    
    lastHeartbeat = undefined;
    
    heartbeat();
  };
  
  // stop sending heartbeat requests
  // and don't emit any "death" events 
  multiplex.die = function() {
    opts.heartbeat = undefined;
  };
  
  // stop responding to heartbeats
  multiplex.playDead = function(state) {
    playingDead = (state === undefined) ? true : false;
  };
  
  return multiplex;
}


rpcMultiStream.syncStream = function(fn, sopts, type) {
  sopts = sopts || {};
  type = type || 'd';
  var i, a;
  if(typeof sopts === 'number') {
    a = [];
    for(i=0; i < sopts; i++) {
      a.push({type: type});
    }
    sopts = a;
  }
  if(!(sopts instanceof Array)) {
    sopts = [sopts];
  }
  for(i=0; i < sopts.length; i++) {
    sopts[i].type = sopts[i].type || type;
    if(sopts[i].type === 'read') {
      sopts[i].type = 'r';
    } else if(sopts[i].type === 'write') {
      sopts[i].type = 'w';
    } else if(sopts[i].type === 'duplex') {
      sopts[i].type = 'd';
    }
  }
  Object.defineProperty(
    fn, 
    '_rpcOpts', 
    {enumerable: false, value: sopts}
  );
  
  return fn;
};


rpcMultiStream.syncReadStream = function(fn, opts) {
  opts = opts || {};
  return rpcMultiStream.syncStream(fn, opts, 'r');
};

rpcMultiStream.syncWriteStream = function(fn, opts) {
  opts = opts || {};
  return rpcMultiStream.syncStream(fn, opts, 'w');
};


module.exports = rpcMultiStream;



