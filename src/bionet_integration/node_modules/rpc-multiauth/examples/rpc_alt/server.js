#!/usr/bin/env nodejs 

var http = require('http');
var path = require('path');

// server static files
var ecstatic = require('ecstatic')({
    root: path.join(__dirname)
});

var websocket = require('websocket-stream');
var rpc = require('rpc-multistream');
var auth = require('../../index.js');

var server = http.createServer(function (req, res) {
    return ecstatic(req, res);
});

var settings = {
  port: 3000,
  host: 'localhost',
  secret: 'my_unique_secret_string' // the token secret
}

server.listen(settings.port, settings.host);

websocket.createServer({server: server}, function(stream) {

    var rpcServer = rpc(auth({ 
        secret: settings.secret,
        tokenExpiration: 4, // expiration in hours
        userDataAsFirstArgument: true, // pass userData as first argument of all rpc functions
        // the login function
        login: function(loginData, cb) {
            if(loginData.username != 'cookiecat' && loginData.password != 'sneeple') {
                return cb("Invalid username or passsword");
            }
            
            // return the user's information
            // second argument must be a unique identifier for the user
            cb(null, 'cookiecat', {
                name: 'cookiecat',
                group: 'user'
            });
        }
        
    }, { // RPC functions
        
        // function with no namespace
        foo: function(userData, cb) {
            cb(null, "hi " + userData.name + " i am foo!");
        },
        
        // functions in the 'user' namespace
        user: {
            bar: function(userData, cb) {
                if(userData.name == 'cookiecat') {
                    cb(null, "I love you Cookie Cat!");
                    return;
                }
                cb(null, "you called the user function mr. " + userData.name);
            }
        }, 
        
        // functions in the 'admin' namespace
        admin: {
            baz: function(userData, cb) {
                cb(null, "you called the admin function");
            }
        }

    // This function is called before all RPC functions.
    // If calls the callback with an error then an authentication
    // error is returned to caller.
    // You can also simply not call the callback and the client
    // will never get a response.
    }, function(userData, namespace, methodName, callback) {
        if(namespace && userData.group != namespace) {
            return callback("Not allowed");
        }
        callback();
    }));

    rpcServer.pipe(stream).pipe(rpcServer);
});

console.log("Server listening on: http://"+settings.host+":"+settings.port+"/");
