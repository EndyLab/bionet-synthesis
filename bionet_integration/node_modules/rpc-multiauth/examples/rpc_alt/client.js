

var rpc = require('rpc-multistream');
var websocket = require('websocket-stream');
var auth = require('../../index.js');


function log(str) {
    var txt = document.getElementById('debug');
    txt.innerHTML += str + "\n";
    txt.scrollTop = txt.scrollHeight;
}

function pageInit() {

    var loginBtn = document.getElementById('login-btn');
    var fooBtn = document.getElementById('foo-btn');
    var barBtn = document.getElementById('bar-btn');
    var bazBtn = document.getElementById('baz-btn');
    var logoutBtn = document.getElementById('logout-btn');
  
    var stream = websocket('ws://' + window.document.location.host);
    var client = rpc();

    client.pipe(stream).pipe(client)

    log("Page loaded");

    client.on('remote', function(remote) {
        log("Connected to server");
        loginBtn.disabled = false;
        fooBtn.disabled = false;
        barBtn.disabled = false;
        bazBtn.disabled = false;
        logoutBtn.disabled = false;
        
        auth.authenticate(remote, function(err, userData) {
            if(err) {
                log("Not logged in");
            } else {
                log("Logged in as: " + userData.name);
            }
        });

        loginBtn.addEventListener('click', function() {
            auth.login(remote, {
                username: 'cookiecat',
                password: 'sneeple'
            }, function(err, token, userData) {
                if(err) return log("Error: " + err);
                log("login successful! token: " + token + " userData: " + JSON.stringify(userData));
            });
        });

        fooBtn.addEventListener('click', function() {
            remote.foo(function(err, msg) {
                if(err) return log("Error: " + err);
                log("foo said: " + msg);
            });
        });

        barBtn.addEventListener('click', function() {
            remote.bar(function(err, msg) {
                if(err) return log("Error: " + err);
                log("bar said: " + msg);
            });
        });

        bazBtn.addEventListener('click', function() {
            remote.baz(function(err, msg) {
                if(err) return log("Error: " + err);
                log("baz said: " + msg);
            });
        });

        logoutBtn.addEventListener('click', function() {
            auth.logout(remote, function(err) {
                if(err) return log("Error: " + err);
                log("Logged out");
            });
            // note: you can also log out when you're disconnected
            // just call logout synchronously and with no options:
            //   auth.logout(); 
            // this deletes auth tokens stored in cookies/localstorage
        });
    });
}


function ready(fn) {
  if (document.readyState != 'loading'){
    fn();
  } else {
    document.addEventListener('DOMContentLoaded', fn);
  }
}

ready(pageInit);
