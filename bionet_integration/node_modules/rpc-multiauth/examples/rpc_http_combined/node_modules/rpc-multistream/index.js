var rpc = require('rpc-stream');
var MuxDemux = require('mux-demux');

function server(methods) {
    var rpcs = rpc(methods, true);
    var mdm = MuxDemux({error: true});

    function genManifest(methods) {
        var manifest = {};
        for(var name in methods) {
            manifest[name] = methods[name]._rpcType || 'async';
        }
        return manifest;
    }

    mdm.on('connection', function(con) {
        con.on('error', function() {});

        if(con.meta == 'rpc') {
            con.pipe(rpcs).pipe(con);
            return;
        }

        if(con.meta == 'meta') {
            con.write(genManifest(methods));
            return;
        }

        try {
            var stream = methods[con.meta[0]](con.meta.slice(1));
            if(stream.readable) {
                stream.pipe(con);
            }
            if(stream.writable) {
                con.pipe(stream);
            }
        } catch(err) {
            mdm.emit('error', err);
        }
    });
    return mdm;
}

function client() {
    var rpcs = rpc(null, true);
    var mdm = MuxDemux({error: true});
    var rpcStream = mdm.createStream('rpc');
    rpcStream.on('error', function() {});
    rpcs.pipe(rpcStream).pipe(rpcs);
    var metaStream = mdm.createStream('meta');

    function createRemoteStream(name, type) {

        return function() {
            var args = [].slice.call(arguments);

            var cb = ('function' == typeof args[args.length - 1])
                ? args.pop()
                : null;
            
            var streamArgs = [name].concat(args);

            var newStream = (
                type === 'readable'
                    ? mdm.createReadStream(streamArgs)
                    : type == 'writable'
                    ? mdm.createWriteStream(streamArgs)
                    : type == 'duplex'
                    ? mdm.createStream(streamArgs)
                    : (function () { throw new Error("unknown stream type: " + type) })()
            )
            newStream.autoDestroy = false;

            return newStream;
        };
    }

    metaStream.on('data', function(manifest) {
        var methods = {};
        for(var name in manifest) {
            if(manifest[name] == 'async') {
                methods[name] = rpcs.createRemoteCall(name);
            } else {
                methods[name] = createRemoteStream(name, manifest[name]);
            }
        }
        mdm.emit('remote', methods);
    }.bind(mdm));
    return mdm;
}


module.exports = function(obj) {
    if(obj) {
        return server(obj);
    }
    return client();
};

module.exports.readable = function(fn, args) {
    fn._rpcType = 'readable';
    return fn;
};

module.exports.writable = function(fn, args) {
    fn._rpcType = 'writable';
    return fn;
};

module.exports.duplex = function(fn, args) {
    fn._rpcType = 'duplex';
    return fn;
};

