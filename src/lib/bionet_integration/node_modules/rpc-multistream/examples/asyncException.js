#!/usr/bin/env node

// This is an example of remote exceptions triggering 
// a local callback call with an error argument

var rpc = require('../index.js');

var server = rpc({

    foo: function(cb) {

        throw new Error("DON'T PANIC");

    }

});

var client = rpc();

client.pipe(server).pipe(client);

client.on('methods', function(methods) {

    methods.foo(function(err) {

        if(err) {
            console.log("Remote error:", err.message);
            return;
        }

        console.log("No error I guess?");
    });
});
