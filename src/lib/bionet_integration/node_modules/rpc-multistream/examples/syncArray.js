#!/usr/bin/env nodejs

// Example of remote functions that return multiple streams as an array

var fs = require('fs');
var rpc = require('../index.js');

var server = rpc({

    multiRead: rpc.syncReadStream(function() {
        var s1 = fs.createReadStream('foo.txt', {encoding: 'utf8'});
        var s2 = fs.createReadStream('foo2.txt', {encoding: 'utf8'});
        return [s1, s2];
        
    // this number indicates how many streams will be returned
    }, 2),


    multiRW: rpc.syncStream(function() {
        var s1 = fs.createReadStream('foo.txt', {encoding: 'utf8'});
        var s2 = fs.createWriteStream('/tmp/foo', {encoding: 'utf8'});
        return [s1, s2];
        
    // differnet ops for each stream
    }, [{
            type: 'read',
            encoding: 'utf8',
            objectMode: false
        },
        {
            type: 'write',
            encoding: 'utf8',
            objectMode: false
        }]),

});

var client = rpc();

client.pipe(server).pipe(client);

client.on('methods', function(methods) {
    var streams = methods.multiRead();
    streams[0].pipe(process.stdout);
    streams[1].pipe(process.stderr);


    var rwStreams = methods.multiRW();
    rwStreams[0].pipe(process.stdout);
    rwStreams[1].write("I am the contents of /tmp/foo\n");

});
