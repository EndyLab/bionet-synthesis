#!/usr/bin/env node

// Example showing asynchronous remote functions that return streams

var path = require('path');
var fs = require('fs');
var through = require('through2');
var rpc = require('../index.js');

var server = rpc({

    // write stream example
    foo: function(filename, cb) {
        var w = fs.createWriteStream(path.join('/tmp', filename), {encoding: 'utf8'});
        cb(null, w);
    },

    // read stream example
    bar: function(cb) {
        var r = fs.createReadStream('foo.txt', {encoding: 'utf8'});
        cb(null, "bar says hi", r);
    }, 

    // multiple streams for one callback
    baz: function(filename, cb) {
        var w = fs.createWriteStream(path.join('/tmp', filename), {encoding: 'utf8'});
        var r = fs.createReadStream('foo.txt', {encoding: 'utf8'});

        cb(null, r, w, {chiao: "hiya!"});
    },

    // duplex stream example
    duper: function(cb) {
        // create duplex transform stream with utf8 input and output

        var ds = through({encoding: 'utf8', decodeStrings: false}, function(data) {
            data = data.toUpperCase();
            this.push(data);
        });

        cb(null, ds);

    }
});

var client = rpc();

client.pipe(server).pipe(client);

client.on('methods', function(methods) {

    methods.bar(function(err, msg, r) {
        console.log("client got: " + msg);
        r.pipe(process.stdout);

        methods.foo('haha.txt', function(err, w) {
            w.write("woop!");
            w.end();

            methods.baz('cookie-cat.txt', function(err, r, w, msg) {
                console.log("remote said:", msg);
                r.on('data', function(data) {
                    console.log("ReadStream: " + data);
                });
                w.write("a treat for your tummy!");
                w.end();

                methods.duper(function(err, dupStream) {


                    console.log("Returned");
                    dupStream.pipe(process.stdout);
                    dupStream.write("making this uppercaaaaase\n");

                });
            });
        });
    });
});
