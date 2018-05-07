var isNode = require('detect-node'); // detect if this is node (or browser)
var jwt = require('jsonwebtoken'); // json web tokens (jwt)
var cookie = require('cookie-cutter'); // cookie parsing and setting
var cookieBaker = require('cookie').serialize; // generate cookie headers
var extend = require('extend');
var clone = require('clone');

/*
  TODO set cookies on server side on login

*/

var defaultTokenName = 'AuthToken';
var defaultTokenHeader = 'authorization';
var defaultTokenExpiration = 336; // in hours (two weeks default)

function unixEpochTime(hours) {
    hours = hours || 0;
    return Math.floor(new Date((new Date()).getTime() + 1000 * 60 * 60 * hours).getTime() / 1000);
}

function createToken(uniqueID, secret, userData, opts) {
    opts = extend({
        expiration: defaultTokenExpiration
    }, opts || {});

    var o = {
        id: uniqueID,
        created: unixEpochTime(),
        expires: unixEpochTime(opts.expiration)
    };

    if(userData) {
        o.userData = userData;
    }

    return jwt.sign(o, secret);
}

function verifyToken(token, secret, opts, cb) {
    opts = extend({
        expiration: defaultTokenExpiration, // how long from creation does token expire (hours)
    }, opts || {});

    jwt.verify(token, secret, function(err, decoded) {
        if(err) return cb("Invalid token: " + err);

        if(decoded.created < unixEpochTime(-opts.expiration)) {
            return cb("Token expired.");
        }

        return cb(null, decoded.userData);
    });
}

function authHTTP(opts) {

    opts = extend(true, {
        allowCookieToken: defaultTokenName, // name of the cookie or false
        allowHeaderToken: defaultTokenHeader, // name of header or false
        tokenExpiration: defaultTokenExpiration, // how long from creation does token expire (hours)
        cookie: {
            setCookie: true, // set cookie on server side
            httpOnly: false, // use httpOnly cookies
            secure: false, // use secure cookies (https only)
            // additional opts:
            //   domain, path, firstPartyOnly, maxAge
        },
        secret: null, // must be set to your unique secret
        success: function(req, res, match, tokenData) {
            if(match && match.next && (typeof match.next == 'function')) {
                match = match.next();
                if(!match) {
                    res.statusCode = 404;
                    res.setHeader("Content-Type", "text/plain");
                    res.end("No matching routes");
                    return;
                }
                match.fn(req, res, match, tokenData);
            } else if(match && (typeof match == 'function')) {
                return match(null, tokenData);
            } else {
                throw new Error(module.name + ": Authentication succeeded but no idea what to do now. Third argument should be a function like next() or an object with a .next() function.");
            }
        },
        fail: function(err, req, res, match) {
            if(match && (typeof match == 'function')) {
                return match(err);
            }
            res.statusCode = 401;
            res.setHeader("Content-Type", "text/plain");
            res.end(err);
        },
        // An optional extra function to run before deciding 
        // if auth check was successful
        // receives args: req, res, match, tokenData, callback
        check: undefined
    }, opts || {});

    if(!opts.secret) {
        throw new Error(module.name + ": You must supply a secret!");
    }
        
    function authFail(errs, req, res, match, cb) {

        var err = errs.join("\n");
        if(!err) {
            if(!opts.allowCookieToken && !opts.allowHeaderToken) {
                err = "Neither cookie nor custom-header authentication allowed so no way to authenticate.";
            } else {
                err = "Unknown error.";
            }
        }

        return cb("Authentication failed: " + err, req, res, match);
    }

    function loginFunc(res, uniqueID, userData, callback) {
        if(typeof res != 'object') {
            callback = userData;
            userData = uniqueID;
            uniqueID = res;
            res = undefined;
        }
        if(typeof userData == 'function') {
            callback = userData;
            userData = null;
        }
        if(!uniqueID || !(typeof uniqueID == 'string' || typeof uniqueID == 'number')) {
            return callback("Missing or invalid ID");
        }   
           
        var token = createToken(uniqueID, opts.secret, userData, {
            expiration: opts.tokenExpiration
        });
        
        if(!opts.allowCookieToken || !res) {
            return callback(null, token);

        }

        var cookieOpts = clone(opts.cookie);
        if(cookieOpts.expiration) {
            cookieOpts.expires = new Date((new Date).getTime() + 1000 * 60 * 60 * opts.tokenExpiration);
        }

        if(res && opts.cookie) {
            res.setHeader('Set-Cookie', cookieBaker(opts.allowCookieToken, token, cookieOpts));
        }

        callback(null, token);
    }

    var authFunc = function(req, res, match) {

        if(typeof res == 'function') {
            match = res;
            res = undefined;
        }

        var successCallback = opts.success;
        var failCallback = opts.fail;

        var token = false;
        var errs = [];

        if(!opts.secret) {
            errs.push("Token secret missing. Cannot authenticate.");
            return authFail(errs, req, res, match, failCallback);
        }
        
        if(opts.allowCookieToken) {
            if(!req.headers.cookie) {
                errs.push("Request had no token cookie.");
            } else {
                var cookies = cookie(req.headers.cookie);
                token = cookies.get(opts.allowCookieToken);
                if(!token) {
                    errs.push("Request did not have the '"+opts.allowCookieToken+"' cookie set.");
                }
            }
        }

        if(!token && opts.allowHeaderToken) {
            if(!req.headers[opts.allowHeaderToken]) {
                errs.push("Request had no token header.");
            } else {
                token = req.headers[opts.allowHeaderToken];
                if(!token) {
                    errs.push("Request did not have the '"+opts.allowHeaderToken+"' custom header set.");
                }
            }
        }

        if(!token) {
            return authFail(errs, req, res, match, failCallback);
        }

        verifyToken(token, opts.secret, {expiration: opts.tokenExpiration}, function(err, tokenData) {

            if(err) {
                errs.push(err);
                authFail(errs, req, res, match, failCallback);
                return;
            }

            if(opts.check && typeof opts.check == 'function') {
                opts.check(req, res, match, tokenData, function(err) {
                    if(err) {
                        errs.push(err);
                        authFail(errs, req, res, match, failCallback);
                        return;
                    }
                    successCallback(req, res, match, tokenData);
                });
            } else {
                successCallback(req, res, match, tokenData);
            }
        });
    };

    var inheritFunc = function(newOpts) {
        var oldOpts = clone(opts);
        newOpts = extend(true, oldOpts, newOpts);
        return authHTTP(newOpts);
    };

    authFunc.login = loginFunc;
    authFunc.inherit = inheritFunc;

    return authFunc;
}


function authRPC(opts, procs, hookOrNamespace) {
    opts = extend({
        tokenExpiration: defaultTokenExpiration,
        userDataAsFirstArgument: false // if true all functions will receive userData as first argument (undefined if user isn't logged in)
    }, opts || {});

    if(!opts.secret) throw new Error("Token auth needs a .secret set");
    if(!opts.login || typeof opts.login != 'function') throw new Error("Token auth needs a .secret set");

    var rpcMethods = {};

    // remember (on the server side) that user is logged
    function rememberUser(token) {
        
        var decoded = jwt.decode(token);

        if(!rpcMethods._rpcMultiAuthData) {

            Object.defineProperty(
                rpcMethods, 
                '_rpcMultiAuthData', 
                {enumerable: false, configurable: true, value: decoded}
            );
        } else {
            rpcMethods._rpcMultiAuthData = decoded;
        }
    }

    // forget (on the server side) user
    function forgetUser() {
        delete rpcMethods._rpcMultiAuthData;
    }

    rpcMethods.login = function(loginData, cb) {

        opts.login(loginData, function(err, id, userData) {
            if(err) return cb("Login failed: " + err);
            var token = createToken(id, opts.secret, userData, {
                expiration: opts.tokenExpiration
            });

            rememberUser(token);

            cb(null, token, userData);
        });
    };

    rpcMethods.logout = function(cb) {
        forgetUser();
        cb();
    };

    // authenticate if client already has token
    rpcMethods.authenticate = function(token, cb) {
        verifyToken(token, opts.secret, {expiration: opts.tokenExpiration}, function(err, decoded) {
            if(err) return cb(err);

            rememberUser(token);
            cb(null, decoded);
        });
    };

    var reservedMethods = Object.keys(rpcMethods);

    function isReserved(methodName) {
        if(!methodName) return false;
        if(reservedMethods.indexOf(methodName) >= 0 || methodName[0] == '_') {
            console.error("rpc method", methodName, "ignored (reserved name)");
            return true;
        }
        return false;
    }

    function authWrap(method, methodName, namespace, hookOrNamespaceMatch) {
        var hook;
        var namespaceMatch;
        
        if(typeof hookOrNamespace == 'string') {
            namespaceMatch = hookOrNamespaceMatch;
        } else if(typeof hookOrNamespaceMatch == 'function') {
            hook = hookOrNamespaceMatch;
        }
        
        function rpcFail(args, err) {
            // if last argument is a function assume 

            if(args.length > 0 && typeof args[args.length-1] == 'function') {
                return args[args.length-1](err);
            } else {
                console.error("Auth failed for synchronous rpc function. Returning null to client but should be throwing exception. TODO implement client side exceptions.");
                return null; // TODO figure out how to throw an exception on the remote side
            }
        }

        function authFail(args, msg) {
            return rpcFail(args, "Unauthorized: " + (msg || "Are you sure you are logged in?"));
        }

        var f = function() {

            var userData;
            if(rpcMethods && rpcMethods._rpcMultiAuthData) {
                userData = rpcMethods._rpcMultiAuthData.userData;
            }

            if(opts.userDataAsFirstArgument) {
                [].splice.call(arguments, 0, 0, userData);
            }
            

            // if this is a namespaced function and namespace matching is on
            if(namespace && namespaceMatch) {

                if(!userData) {
                    return authFail(arguments, "Not logged in.");
                }

                if(userData[namespaceMatch] != namespace) {
                    return authFail(arguments, "You are not allowed to call procedures in the '" + namespace + "' namespace");
                }

                return method.apply(this, arguments);
                
            } else {
                // if user specified a hook function instead of namespace matching
                if(hook) {
                    var args = arguments;
                    hook(userData, namespace, methodName, function(err) {
                        if(err) return authFail(args, err);
                        return method.apply(this, args);
                    });
                } else {
                    return method.apply(this, arguments);
                }
            }
        };

        // Propagate synchronous function tag.
        // This tag ensures that the function is recognized
        // as the appropriate read/write/duplex stream returning function
        if(method._rpcOpts) {
            Object.defineProperty(
                f, 
                '_rpcOpts', 
                {enumerable: false, value: method._rpcOpts}
            );
        }

        return f;
    }

    var key, innerKey, namespace, wrapped;
    for(key in procs) {
        if(isReserved(key)) continue;

        if(typeof procs[key] == 'function') {

            rpcMethods[key] = authWrap(procs[key], key, null, hookOrNamespace);

        } else if(typeof procs[key] == 'object') {
            namespace = key;
            for(innerKey in procs[key]) {
                if(typeof procs[key][innerKey] != 'function') continue;
                if(isReserved(innerKey)) continue;

                rpcMethods[innerKey] = authWrap(procs[key][innerKey], innerKey, namespace, hookOrNamespace);

            }
        }
    }
    
    return rpcMethods;
}

var rpcMultiAuth;


var serverExport = rpcMultiAuth = function() {

        // Being initialized for RPC auth
        if(typeof arguments[0] == 'object' && typeof arguments[1] == 'object') {
            var opts = arguments[0];
            var procs = arguments[1]; // remote procedures
            var hookOrNamespace = arguments[2];
            
            return authRPC(opts, procs, hookOrNamespace);

        // Being initialized for HTTP auth (e.g. as middleware)
        } else if((typeof arguments[0] == 'string' || typeof arguments[0] == 'object')) {
            var opts;
            var callback;
            if(typeof arguments[0] == 'object') {
                opts = arguments[0];
            } else {
                opts = {
                    secret: arguments[0]
                };
            }
            callback = arguments[1];

            return authHTTP(opts, callback);
            
        } else {
            throw new Exception("Wrong arguments. Read the docs.");
        }
    };



    var store = require('store'); // LocalStorage wrapper
    var xhr = require('xhr'); // XMLHttpRequest wrapper
    
var clientExport = rpcMultiAuth = {

        // save token using cookie if opts.setCookie is true (default false)
        // and save in LocalStorage if opts.useLocalStorage is true (default true)
        saveToken: function(token, opts) {

            opts = extend({
                tokenName: defaultTokenName, // cookie name and/or localstorage key to use
                setCookie: false, // whether to save token in cookie on client side
                useLocalStorage: true // true to save token using LocalStorage
            }, opts || {});
            var saved = false;

            if(opts.setCookie) {
                var decoded = jwt.decode(token);
                if(decoded && decoded.expires) {
                    document.cookie = cookieBaker(opts.tokenName, token, {
                        expires: new Date(decoded.expires * 1000)
                    });
                    saved = true;
                }
            }
            if(opts.useLocalStorage) {
                store.set(opts.tokenName, token);
                saved = true;
            }
            return saved;
        },

        // do a purely local check to see if we have a token
        // and if it is still good (not expired)
        isLoggedIn: function(opts) {
            var token = rpcMultiAuth.getToken(opts);

            if(!token) return false

            var decoded = jwt.decode(token);

            if(!decoded || !decoded.expires) {
                return false;
            }
            if(decoded.expires > unixEpochTime()) {
                return true;
            }
            return false;
        },

        // delete the token
        delToken: function(opts) {
            opts = extend({
                tokenName: defaultTokenName, // cookie name and/or localstorage key to use
                delCookie: true, // whether to delete token in cookie on client side
                useLocalStorage: true // true to delete token using LocalStorage
            }, opts || {});

            var removed = false;
            if(opts.delCookie) {
                document.cookie = cookieBaker(opts.tokenName, '', {
                    expires: new Date(0)
                });
                removed = true;
            }
            if(opts.useLocalStorage) {
                store.remove(opts.tokenName);
                removed = true;
            }
            return removed;
        },

        logout: function(opts, callback) {
            var remote;
            
            // if this is called with an rpc object that supports logouts
            if(typeof opts == 'object' && typeof opts.logout == 'function') {
                remote = opts;
                rpcMultiAuth.delToken();
                remote.logout(callback);
            } else {
                return rpcMultiAuth.delToken(opts);
            }
        },

        // get token from cookie if available 
        // otherwise from localstorage if available
        getToken: function(opts) {

            opts = extend({
                tokenName: defaultTokenName, // cookie name and/or localstorage key to use
                token: undefined // of opts contains a .token then that will be used
            }, opts || {});

            if(opts.token) return opts.token;

            var token;
            var cookies = cookie(document);

            if(!token) {
                token = cookie.get(opts.tokenName);
            }
            if(!token) {
                token = store.get(opts.tokenName);
            }
            return token;
        },
        
        // check if token is valid / if we are authorized
        // remote can be an rpc-multiauth (or other rpc) instance 
        // or a URL to use for HTTP POST login
        // if it is null then the relative url 'login' will be used
        authenticate: function(remote, opts, callback) {
            if(typeof opts == 'function') {
                callback = opts;
                opts = null;
            }
            opts = opts || {};
            if(typeof opts == 'string') {
                opts = {
                    token: opts
                };
            }
            opts = extend({
                token: undefined, // optionally explicitly pass a token
                setCookie: false, // whether to save token in cookie on client side
                tokenName: defaultTokenName, // cookie name and/or localstorage key to use
                useLocalStorage: true, // true to save token using LocalStorage
                httpMethod: 'POST', // which method to use if using http requests
                tokenHeader: defaultTokenHeader // header to use for token, set to false to disable sending token in header (e.g. if using cookies)
            }, opts);

            remote = remote || 'authenticate'; // default URL to use for auth
                
            var token = rpcMultiAuth.getToken(opts);
            if(!token) {
                return callback('No token available so could not authenticate. Either explicitly pass a token in the opts argument like so: {token: "my_token"} or ensure that you have a token saved using cookies or LocalAuth using the cookie name / localstorage key "'+opts.tokenName+'".');
            }

            if(typeof remote == 'string') {
                var o = {
                    uri: remote,
                    method: opts.httpMethod.toUpperCase(),
                    headers: {}
                };

                if(opts.tokenHeader) {
                   o.headers[opts.tokenHeader] = token;
                }
                xhr(o, function(err, resp, body) {
                        if(err || resp.statusCode < 200 || resp.statusCode >= 300) {
                            return callback(err || body);
                        }
                    if(!body) return callback("Got empty response from server but no http error status code.");
                    if(typeof body != 'string') return callback("Got non-string response from server. Expected userData.");
                    callback(null, body);
                });
            } else {
                remote.authenticate(token, callback);
            }
        },
        
        // log in
        // remote can be an rpc-multiauth (or other rpc) instance 
        // or a URL to use for HTTP POST login
        // if it is null then the relative url 'login' will be used
        login: function(remote, loginData, opts, callback) {
            if(typeof opts == 'function') {
                callback = opts;
                opts = null;
            }
            opts = extend({
                setCookie: false, // whether to save token in cookie on client side
                useLocalStorage: true, // whether to save token in LocalStorage
                tokenName: defaultTokenName, // cookie name and/or localstorage key to use
                tokenHeader: defaultTokenHeader // header to use for token, set to false to disable sending token in header (e.g. if using cookies)
            }, opts || {});

            remote = remote || 'login'; // default URL to use for login

            function doLogin(remote, loginData, opts, callback) {
                if(typeof remote == 'string') {
                    xhr({
                        uri: remote,
                        method: 'POST',
                        body: JSON.stringify(loginData)
                    }, function(err, resp, token) {
                        if(err || resp.statusCode < 200 || resp.statusCode >= 300) {
                            return callback(err || token);
                        }
                        if(!token) return callback("Got empty response from server but no http error status code.");
                        if(typeof token != 'string') return callback("Got non-string response from server. Expected token.");
                        
                        rpcMultiAuth.saveToken(token, opts);
                        
                        callback(null, token);
                    });
                } else {
                    remote.login(loginData, function(err, token, userData) {
                        if(err) return callback(err);
                        rpcMultiAuth.saveToken(token, opts);
                        callback(null, token, userData);
                    });
                }
            }

            doLogin(remote, loginData, opts, callback);

        },

        // return a wrappted xhr requester
        requester: function(popts) {

            if(typeof popts == 'string') {
                popts = {
                    token: popts
                }
            }

            popts = extend({
                tokenHeader: defaultTokenHeader, // name of header field to use for auth token. if false, disable sending of auth token using its own header field (rely on cookies only)
                tokenName: defaultTokenName, // name of cookie and/or localstorage field to look for
                token: undefined // explicitly set an auth token to use. if not set then requester will look for token in cookie
            }, popts || {});

            return function(opts, callback) {
                opts.headers = opts.headers || {};
                if(popts.tokenHeader) {
                    var token = opts.token || popts.token;
                    token = rpcMultiAuth.getToken({token: token});
                    if(token) {
                        opts.headers[popts.tokenHeader] = token;
                    }
                }
                
                return xhr(opts, callback);
            };
        }

    };

    rpcMultiAuth.xhr = xhr;


// For backwards compatibility with older versions of rpc-multiauth
if(isNode) { // server side
  serverExport.server = serverExport;
  serverExport.client = clientExport;
  module.exports = serverExport;
} else { // client side
  clientExport.server = serverExport;
  clientExport.client = clientExport;
  module.exports = clientExport;
}
