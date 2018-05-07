#!/usr/bin/env nodejs 

var http = require('http');
var path = require('path');

// server static files
var ecstatic = require('ecstatic')({
    root: path.join(__dirname)
});

var websocket = require('websocket-stream');
var rpc = require('rpc-multistream');

// server static files
var ecstatic = require('ecstatic')({
    root: path.join(__dirname)
});

var auth = require('../../index.js');

var settings = {
  port: 3000,
  host: 'localhost',
  secret: 'my_unique_secret_string' // the token secret
}

var myAuth = auth({
  secret: settings.secret,
  cookie: {
      setCookie: false // don't set cookie on server side
  }
});

var server = http.createServer(function (req, res) {

    if(req.url == '/' || req.url == '/bundle.js') {
        return ecstatic(req, res);

    } else if(req.url == '/private') {
        res.setHeader("Content-Type", "text/plain");

        myAuth(req, function(err, tokenData) {
            if(err) {
                res.statusCode = 401;
                res.end("Unauthorized: " + err);
                return;
            }
            res.end("Yay you are authorized!");
        });
    } else {
        myAuth(req, function(err, tokenData) {
            if(err) {
                res.statusCode = 401;
                res.end("Unauthorized: " + err);
                return;
            }
            return ecstatic(req, res);
        });
    }
});

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
        foo: function(cb) {
            cb(null, "hi i am foo!");
        },
        
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
