#!/usr/bin/env nodejs

// Examples of using both binary, string and object streams

var fs = require('fs');
var from = require('from2');
var rpc = require('../index.js');

var server = rpc({
  
  // function returning a text read stream
  foo: rpc.syncReadStream(function() {
    var s = fs.createReadStream('foo.txt');
    return s;
  }, {
    encoding: null,
    objectMode: false
  }),

  // function returning a text read stream
  bar: rpc.syncReadStream(function() {
    var s = fs.createReadStream('foo.txt');
    return s;
  }, {
    encoding: 'utf8',
    objectMode: false
  }),

  // function returning a text read stream
  baz: rpc.syncReadStream(function() {
    var i = 0;
    return from.obj(function(size, next) {
      if(i++) return next(null, null);
      next(null, {
        hoopy: 'frood'
      });
    });
  }, {
    objectMode: true
  })
  
}, {
  objectMode: false
//  debug: true
});

var client = rpc(undefined, {
  objectMode: false,
  explicit: true // with explicit set it doesn't matter if defaults differ
});

client.pipe(server).pipe(client);

client.on('methods', function(methods) {

  var binaryStream = methods.foo();

  binaryStream.on('data', function(data) {
    console.log("binaryStream got", typeof data, ':', data);
  });

  var utf8Stream = methods.bar();

  utf8Stream.on('data', function(data) {

    console.log("utf8Stream got", typeof data, ':', data);

  });

  var objStream = methods.baz();

  objStream.on('data', function(data) {

    console.log("objStream got", typeof data, ':', data);

  });
});
