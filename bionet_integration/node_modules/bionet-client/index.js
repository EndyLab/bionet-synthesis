
var auth = require('rpc-multiauth').client;
var path = require('path');
var fs = require('fs');
var websocket = require('websocket-stream');
var rpc = require('rpc-multistream');
var extend = require('xtend');

module.exports = function(url, opts, callback) {
  if(typeof opts === 'function') {
    callback = opts;
    opts = {}
  }

  url = url.replace(/^https:/, 'wss:');
  url = url.replace(/^http:/, 'ws:');

  opts = extend({
    username: undefined,
    password: undefined
  }, opts || {});

  var stream = websocket(url);

  stream.on('error', function(err) {
    if(exiting) return;
    callback(new Error("connection error"));
  });

  stream.on('close', function() {
    if(exiting) return;
    callback(new Error("connection closed"));
  });

  var rpcClient = rpc(null, {
    objectMode: true,
    heartbeat: 2000
  });

  rpcClient.pipe(stream).pipe(rpcClient);

  rpcClient.on('error', function(err) {
    if(exiting) return;
    callback(err);
  });

  var exiting;
  function done() {
    exiting = true;
    stream.end();
  }

  rpcClient.on('methods', function (remote) {

    if(!opts.username) return callback(null, done, remote);

    auth.login(remote, {
      username: opts.username,
      password: opts.password
    }, function(err, token, userData) {
      if(err) {
        return callback(err);
      }

      return callback(null, done, remote, userData, token);
    });
  });
};
