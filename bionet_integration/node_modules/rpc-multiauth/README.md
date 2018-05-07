
rpc-multiauth is simple token-based authentication and access control for both RPC over http sockets and traditional http queries.

rpc-multiauth provides optional cookie-based authentication, still using tokens, to support authenticated traditional queries when client side control over headers isn't possible (e.g. automatic loading of images for <img> tags)..

rpc-multiauth also provides client-side helper methods for authenticated XMLHTTPRequests without cookies by wrapping the [xhr](https://www.npmjs.com/package/xhr) module.

rpc-multiauth works especially well with the npm modules [routes](https://www.npmjs.com/package/routes) and [rpc-multistream](https://www.npmjs.com/package/rpc-multistream), but can be used alone or with other modules.

# RPC authentication

## Server side

```
var auth = require('rpc-multiauth').server;

var rpcServer = rpc(auth(
  { // the following options _must_ be set
    secret: 'my_unique_secret_string',
    login: function(data, cb) {
      if(!isValidUserAndPassword(data)) return cb("Invalid username or passsword");      cb();
    }

  }, { // RPC functions
    // no namespace
    foo: function(cb) {
      console.log("foo called");
      cb(null, "foo says hi");
    },

    // namespaced functions
    user: {
      bar: function(cb) {

      }
    }, 

    // more namespaced functions
    admin: {
      baz: function(cb) {

      }
    }


  // can be either:
  // * undefined: In which case all namespaced functions just require login
  // * string: In which case user[string] must match the namespace name (which can be a regex)
  // * function: receives (userData, namespace, functionName, cb) and calls back with an error if access is denied or nothing/null/undefined if access is granted
  }, function(userData, namespace, functionName, callback) {

  });
```

When auth fails for async functions, if the last argument is a function, it is assumed to be a callback with the first argument being an optional error and is called with an "Unauthorized" error.

Note: Currently there is no real error reporting when auth fails for a sync function. The server will spit out a console.error and null will be returned. We need to solve [this issue](https://github.com/Juul/rpc-multistream/issues/2) so we can throw exceptions that propagate to the client.

## Client side

Simlpe usage:

```
var auth = require('rpc-multiauth').client;

function rpcConnected(remote) {

  auth.authenticate(remote, function(err, userData) {
    if(err) {
      console.log("Not logged in");
    } else {
      console.log("Logged in as: " + userData.name);
    }
  });

  loginBtn.addEventListener('click', function() {

    auth.login(remote.login, {
      user: 'marc@juul.io',
      password: 'foo'
    }, function(err, token, userData) {
       // do stuff 
    });

  });
}
```

With opts:

```
var auth = require('rpc-multiauth').client;

function rpcConnected(remote) {

  auth.authenticate(remote, {
      // These are the default options:
      // Assuming token is not set
      // auth.authenticate will check LocalStorage and cookies for the auth token
      token: undefined, // optionally explicitly pass a token
      setCookie: false, // whether to save token in cookie on client side
      tokenName: 'AuthToken', // cookie name and/or localstorage key to use
      useLocalStorage: true, // true to save token using LocalStorage
      httpMethod: 'POST', // which method to use if using http requests
      tokenHeader: 'authorization' // header to use for token, set
    }, function(err, userData) {
    if(err) {
      console.log("Not logged in");
    } else {
      console.log("Logged in as: " + userData.name);
    }
  });

  loginBtn.addEventListener('click', function() {

    auth.login(remote.login, {'marc@juul.io', 'foo'}, { 
      // optional opts argument with defaults shown
      setCookie: false, // save the cookie on the client side
      useLocalStorage: true, // whether to save token in LocalStorage
      tokenName: 'Token', // cookie name and/or localstorage key to use
      tokenHeader: 'authorization' // header to use for token, set
    }, function(err, token, userData) {
       // do stuff 
    });

  });
}
```

See examples/ for more info.

# HTTP request authentication

## Server side

Basic usage:

```
var myAuth = auth(settings.secret).server;

http.createServer(function(req, res) {

    if(req.url == '/login') {
        res.setHeader("Content-Type", "text/plain");

        // Your own function to check if the login info is valid
        checkLogin(req, function(err, userData) {
            if(err) {
                res.statusCode = 400;
                res.end("Login error: " + err);
                return;
            }

            // log the user in
            myAuth.login(res, userData.id, userData, function(err, token) {
                if(err) {
                    res.statusCode = 400;
                    res.end("Login error: " + err);
                    return;
                }
                res.end(token);
            });
        })    
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
    }
});
```

Initializing auth for HTTP with options:

```
var myAuth = auth({
  secret: null, // must be set to your unique secret
  tokenExpiration: 336, // how long from creation does token expire (hours)
  cookie: {
    setCookie: true, // set cookie on server side
    httpOnly: false, // use httpOnly cookies
    secure: false, // use secure cookies (https only)
    // additional opts:
    // domain: restrict auth cookies to this domain
    // path: restrict auth cookies to this path
    // firstPartyOnly: see RFC6265
    // maxAge: relative max age of the cookie from when the client receives it
  },
  allowCookieToken: 'AuthToken', // name of the cookie or false
  allowHeaderToken: 'authorization' // name of header or false
});
```

Use as middleware with the routes package:

```
var myAuth = auth('MY_TOKEN_SECRET');

// Assume that login is implemented similarly to previous examples

// Routes 
myrouter.addRoute('/users-only/*?', myAuth);

myrouter.addRoute('/users-only/profile', function(req, res, match, userData) {
  res.end("Profile page for user: " + userData.name);
});
```

Or with multiple auth levels:

```
var userAuth = auth('MY_TOKEN_SECRET');

// inherit options from userAuth but add a check
var adminAuth = userAuth.inherit({
  check: function(req, res, match, userData, callback) {
    if(userData.group != 'admin') {
      return callback("You are not an admin");
    }
    callback();
  }
});

myrouter.addRoute('/', function(req, res, match) {
  res.end("Welcome to the main page!");
});

// Verify that user is logged in
myrouter.addRoute('/users-only/*?', userAuth);
myrouter.addRoute('/users-only/profile', function(req, res, match, userData) {
  res.end("Profile page for user: " + userData.name);
});

// Verify that user is logged in and is an admin
myrouter.addRoute('/admins-only/*?', adminAuth);
myrouter.addRoute('/admins-only/settings', function(req, res, match, userData) {
  res.end("Insert admin settings page here");
});

// Default route
myrouter.addRoute('/*', function(req, res, match) {
  res.statusCode = 404;
  res.end("Page not found");
});

```

## Client side

If you haven't disabled cookies, then after logging in you will simply be able to access protected URLs with normal HTTP requests.

If you have disabled cookies (or simply prefer not to use it), then you can instantiate an http requester that will automatically send the token in a custom authentication header with each request:

```
var req = auth.requester();
```

The simple usage above will automatically find any auth tokens saves as Cookies or in LocalStorage.


or with options:

```
var req = auth.requester({ 
  tokenHeader: 'Authorization', // name of header field to use for auth token. if false, disable sending of auth token using its own header field (rely on cookies only)
  tokenName: 'AuthToken', // name of cookie and/or localstorage field to look for
  token: undefined // explicitly set an auth token to use. if not set then requester will look for token in cookie and LocalStorage
});
```

The result req function is simply a wrapped version of (xhr)[https://www.npmjs.com/package/xhr] that sends along the auth token in a header (unless tokenHeader was set to false) and optionally takes the additional option:

```
  {
    token: 'my_token'
  }
```

in case you want to explicitly define a per-request token.

Example usage: 

```
var auth = require('rpc-multiauth').client;

var req = auth.requester({json: true});

req.post({
  uri: '/intergalactic',
  body: {
    cookie: 'cat'
  }
}, function(err, resp, body) {
   if(err) return console.error("Request failed:", err);
   console.log("Got response:", body);
});
```

# Examples

There are a bunch of complete examples in the `examples/` dir. 

To try an example do s:

```
cd examples/rpc_simple/
npm install
npm run build
npm start
```

Then open `http://localhost:3000/` in your browser.