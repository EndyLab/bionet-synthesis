#!/usr/bin/env nodejs

// Example of using opts.explicit 
// Try turning off opts.explicit on the client and see it fail

var fs = require('fs');
var from = require('from2');
var rpc = require('../index.js');


var server = rpc({
  
  // function returning a text read stream
  foo: rpc.syncReadStream(function() {
    return fs.createReadStream('foo.txt', {encoding: 'utf8'});
  }),
  
  // function returning an object read stream
  bar: rpc.syncReadStream(function() {
    var i = 0;
    var s = from.obj(function(size, next) {
      if(i++) return next(null, null);
      next(null, {
        hoopy: 'frood'
      });
    });
    s.on('error', function(err) {
      console.log(err);
    });
  }, {
    objectMode: true // explicitly specify that this stream is objectMode
  })
}, {
  objectMode: false // default for streams is non-objectMode
});

var client = rpc(undefined, {
  objectMode: true, // different default than the server
  explicit: true // with explicit set it doesn't matter if defaults differ
});

client.pipe(server).pipe(client);

client.on('methods', function(methods) {

  var stringStream = methods.foo();
  stringStream.pipe(process.stdout);
  
  stringStream.on('error', function(err) {
    console.log("stringStream error:", err);
  });

  var objStream = methods.bar();

  objStream.on('data', function(data) {

    console.log("got", typeof data, ':', data);

  });

});
