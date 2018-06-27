
var bc = require('bionet-client');
var test = require('tape')

test('connect', function(t) {

  t.plan(1);

  try {
    var settings = require('../settings.js');
  } catch(e) {
    t.fail("Your settings.js is missing, invalid or unreadable");
    return;
  }

  bc(settings.node_url, function(err, done, remote) {
    if(err) return t.fail("failed to connect: " + err);

    // check that RPC is working by calling the "foo -> bar" function
    remote.foo(function(err, val) {
      if(err) return t.fail("remote.foo call failed: " + err);

      t.equal(val, "bar");

      t.end();
      process.exit(0);
    });
  });
});
