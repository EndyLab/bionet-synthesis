#!/usr/bin/env nodejs

// Example of error reporting for synchronous calling using rpc-multistream
//
// This example is split up into a server and a client since the error reporting
// does not work correctly when both server and client run in the same process

var fs = require('fs');
var net = require('net');
var rpc = require('../index.js');

var server = rpc({

    getStream: rpc.syncStream(function() {
        throw new Error("Something bad happened");
    }),

    // this function has no way of reporting errors back to the caller
    // so the error will get emitted on the rpc stream itself
    // on both ends
    bad: function() {
        throw new Error("Something bad happened");
    }

});

server.on('error', function(err, isRemote, functionName) {
    if(isRemote) {
        console.log("Remote function '"+functionName+"' experienced an exception but had no way to report it back to the caller: " + err);
    } else if(isRemote === false){
        console.log("Local function '"+functionName+"' called by remote experienced an exception but had no way to report it to the caller: " + err);
    } else {
        console.log(err);
    }
});

net.createServer(function (con) {
    con.pipe(server).pipe(con);
}).listen(4242);

console.log("Server listening...");
