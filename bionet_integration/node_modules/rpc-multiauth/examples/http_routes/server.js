#!/usr/bin/env nodejs 

var http = require('http');
var path = require('path');

// server static files
var ecstatic = require('ecstatic')({
    root: path.join(__dirname)
});

var auth = require('../../index.js');
var Routes = require('routes');

var settings = {
  port: 3000,
  host: 'localhost',
  secret: 'my_unique_secret_string' // the token secret
}

function checkLogin(req, callback) {
    var body = '';

    req.on('data', function(data) {
        body += data;
    });

    req.on('end', function() {
        var o = JSON.parse(body);
        if(!o || o.username != 'cookie' || o.password != 'cat') {
            return callback("Invalid username or password");
        }
        callback(null, {
            id: '1234',
            name: 'cookiecat',
            desc: 'A treat for your tummy!',
            group: 'user'
        })
    });
}

var userAuth = auth(settings.secret);

// inherit options from userAuth but add a check
var adminAuth = userAuth.inherit({
  check: function(req, res, match, userData, callback) {
    if(userData.group != 'admin') {
      return callback("You are not an admin");
    }
    callback();
  }
});

var router = new Routes();

router.addRoute('/', ecstatic);
router.addRoute('/bundle.js', ecstatic);

router.addRoute('/public', function(req, res, match) {
    res.setHeader("Content-Type", "text/plain");
    res.end("Public information");
});

router.addRoute('/login', function(req, res, match) {
    res.setHeader("Content-Type", "text/plain");
    console.log("got to login");

    checkLogin(req, function(err, userData) {
        if(err) {
            res.statusCode = 400;
            res.end("Login error: " + err);
            return;
        }
        userAuth.login(res, userData.id, userData, function(err, token) {
            if(err) {
                res.statusCode = 400;
                res.end("Login error: " + err);
                return;
            }
            res.end(token);
        });
    });
});

// Verify that user is logged in
router.addRoute('/user/*?', userAuth);
router.addRoute('/user/profile', function(req, res, match, userData) {
    res.setHeader("Content-Type", "text/plain");
    res.end("Profile info for user: " + userData.name);
});

// Verify that user is logged in and is an admin
router.addRoute('/admin/*?', adminAuth);
router.addRoute('/admin/settings', function(req, res, match, userData) {
    res.setHeader("Content-Type", "text/plain");
    res.end("Admin settings");
});

// Default route
router.addRoute('/*', function(req, res, match) {
    res.statusCode = 404;
    res.setHeader("Content-Type", "text/plain");
    res.end("Page not found");  
});

var server = http.createServer(function(req, res) {
  var match = router.match(req.url);
  match.fn(req, res, match);
});

server.listen(settings.port, settings.host);

console.log("Server listening on: http://"+settings.host+":"+settings.port+"/");
