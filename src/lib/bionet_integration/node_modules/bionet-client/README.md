
This is a client library for connecting to [a bionet node](https://github.com/biobricks/bionet-new).

# Usage

To establish an unauthenticated RPC connection:

```
var bc = require('bionet-client');

bc("https://endylab.stanford.edu", function(err, done, remote) {
  if(err) return console.error(err);

  remote.doSomething(function(err, result) {
    if(err) return console.error(err);

    console.log("the server said:", result);

    done(); // close RPC connection
  });

});
```

To establish an authenticated RPC connection simply supply username and password:

```
bc("https://endylab.stanford.edu", {
  username: 'myuser',
  password: 'mypassword'
}, function(err, done, remote, userData, token) {

  // note that we're now receving the extra arguments:
  // * userData: object containing info about logged-in user
  // * token: authentication token

});
```

The connection is established using [rpc-multistream](https://github.com/biobricks/rpc-multistream) which allows calling remote synchronous functions that return streams and functions with asynchronous callbacks that give streams as arguments.

The authentication is handled by [rpc-multiauth](https://github.com/biobricks/rpc-multiauth).

# Testing

Beware that testing requires a running bionet node to test against.

# Copyright and license

Copyright 2016-2018 BioBricks Foundation

License: AGPLv3