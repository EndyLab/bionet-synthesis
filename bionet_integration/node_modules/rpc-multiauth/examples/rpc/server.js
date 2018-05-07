#!/usr/bin/env nodejs 

var http = require('http');
var path = require('path');
var fs = require('fs');

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
        foo: function(msg, cb) {
            cb(null, 'oh hello! thanks for saying "'+msg+'"');
        },

        myStream: rpc.syncReadStream(function() {
            return fs.createReadStream('myfile.txt', {encoding: 'utf8'});
        }),
        
        // functions in the 'user' namespace
        user: {
            bar: function(cb) {
                cb(null, "you called the user function");
            }
        }, 
        
        // functions in the 'admin' namespace
        admin: {
            baz: function(cb) {
                cb(null, "you called the admin function");
            }
        }

    // match the returned user data's .group to the allowed namespaces
    }, 'group'))

    rpcServer.pipe(stream).pipe(rpcServer);
});

console.log("Server listening on: http://"+settings.host+":"+settings.port+"/");
