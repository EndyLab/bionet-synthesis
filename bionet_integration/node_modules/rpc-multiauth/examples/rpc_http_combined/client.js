
var rpc = require('rpc-multistream');
var websocket = require('websocket-stream');
var auth = require('../../index.js');

var req = function() {log("Not logged in");};

function log(str) {
    var txt = document.getElementById('debug');
    txt.innerHTML += str + "\n";
    txt.scrollTop = txt.scrollHeight;
}

function pageInit() {

    var loginBtn = document.getElementById('login-btn');
    var xhrBtn = document.getElementById('xhr-btn');
    var imageBtn = document.getElementById('image-btn');
    var logoutBtn = document.getElementById('logout-btn');
  
    var stream = websocket('ws://' + window.document.location.host);
    var client = rpc();

    client.pipe(stream).pipe(client)

    log("Page loaded");

    client.on('remote', function(remote) {
        log("Connected to server");
        loginBtn.disabled = false;
        xhrBtn.disabled = false;
        imageBtn.disabled = false;
        logoutBtn.disabled = false;
        
        auth.authenticate(remote, {
            setCookie: true // set cookie on client side
        }, function(err, userData) {
            if(err) {
                log("Not logged in");
            } else {
                log("Logged in as: " + userData.name);
                req = auth.requester();
            }
        });

        loginBtn.addEventListener('click', function() {
            auth.login(remote, {
                username: 'cookiecat',
                password: 'sneeple'
            }, {
                setCookie: true // set cookie on client side
            }, function(err, token, userData) {
                if(err) return log("Error: " + err);

                req = auth.requester();
                log("login successful! token: " + token + " userData: " + JSON.stringify(userData));
            });
        });

        xhrBtn.addEventListener('click', function() {
            req({
                uri: '/private'
            }, function(err, res, body) {
                if(err) return log("Error: " + err + ' | ' + body);
                log("Server said: " + body);
            });
        });

        imageBtn.addEventListener('click', function() {
            var img = document.createElement("IMG");
            img.src = "/example.jpg";
            log("loading image: " + img.src);
            img.onload = function(e) {
                document.body.appendChild(img);
                log("image loaded");
            }
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
