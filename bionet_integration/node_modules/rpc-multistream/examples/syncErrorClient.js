#!/usr/bin/env nodejs

// Example of error reporting for synchronous calling using rpc-multistream
//
// This example is split up into a server and a client since the error reporting
// does not work correctly when both server and client run in the same process

var net = require('net');
var rpc = require('../index.js');

var client = rpc();

client.on('methods', function(methods) {

    // this function has no callback 
    // so it has to report remote errors 
    // by emitting an error on the stream
    var stream = methods.getStream();

    stream.on('error', function(err, isRemote, functionName) {
        if(isRemote) {
            console.log("The remote function '"+functionName+"' experienced an error: " + err);
        } else {
            console.log(err);
        }
    });

    // this function has no way of reporting errors 
    // other than to emit a general error on the 'client' stream
    // (see below)

});

client.on('error', function(err, isRemote, functionName) {
    if(isRemote) {
        console.log("Remote function '"+functionName+"' experienced an exception but had no way to report it back to the caller: " + err);
    } else if(isRemote === false){
        console.log("Local function '"+functionName+"' called by remote experienced an exception but had no way to report it to the caller: " + err);
    } else {
        console.log("Error: " + err);
    }
});

console.log("Client connecting...");
var con = net.connect({port: 4242}, function() {
    console.log("Client connected!");
    con.pipe(client).pipe(con);

});
