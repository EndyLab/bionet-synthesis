
var auth = require('../../index.js');
var xhr = require('xhr'); // XMLHTTPRequests
var jwt = require('jsonwebtoken'); // json web tokens (jwt)


function log(str) {
    var txt = document.getElementById('debug');
    txt.innerHTML += str + "\n";
    txt.scrollTop = txt.scrollHeight;
}

var req = function() {log("Not logged in.")};


function pageInit() {

    var loginBtn = document.getElementById('login-btn');
    var publicBtn = document.getElementById('public-btn');
    var userBtn = document.getElementById('user-btn');
    var logoutBtn = document.getElementById('logout-btn');

    if(auth.isLoggedIn()) {
        log("Already logged in!");
        req = auth.requester();
    } else {
        log("Not logged in");
    }

    loginBtn.addEventListener('click', function() {
        console.log("Clicked the login button");
        auth.login('/login', {
            username: 'cookie',
            password: 'cat'
        }, {
            setCookie: true, // set cookies from client side
            useLocalStorage: false // don't save auth token using LocalStorage
        }, function(err, token) {
            if(err) return log("Error: " + err);
            log("Logged in!");

            req = auth.requester({
                tokenHeader: false // don't send token in custom header (cookie only)
            });

            // If you didn't want to use localstorage nor cookies
            // then you could also manually pass the token like so:
            // req = auth.requester(token);
        });
    });

    publicBtn.addEventListener('click', function() {
        xhr({
            uri: '/public'
        }, function(err, res, body) {
            if(err) return log("Error: " + err + ' | ' + body);
            log("Server said: " + body);
        });
    });
    
    userBtn.addEventListener('click', function() {
        req({
            uri: '/private'
        }, function(err, res, body) {
            if(err) return log("Error: " + err + ' | ' + body);
            log("Server said: " + body);
        });
    });

    logoutBtn.addEventListener('click', function() {
        console.log("Clicked the login button");
        auth.logout();
        log("Logged out!");

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
