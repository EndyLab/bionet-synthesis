
var rpc = require('rpc-multistream');
var websocket = require('websocket-stream');

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
        
        loginBtn.addEventListener('click', function() {
            remote.login({
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
            remote.logout(function(err) {
                if(err) return log("Error: " + err);
                log("Logged out");
            });
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
