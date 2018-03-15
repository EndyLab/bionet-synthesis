
var bc = require('bionet-client');
var test = require('tape')

test('authenticate', function(t) {

  t.plan(1);

  try {
    var settings = require('../settings.js');
  } catch(e) {
    t.fail("Your settings.js is missing, invalid or unreadable");
    return;
  }

  bc(settings.node_url, {
    username: settings.username,
    password: settings.password
  } ,function(err, done, remote) {
    if(err) return t.fail("failed to connect: " + err);

    // check that authenticated RPC is working 
    // by calling the "foo_user -> bar_user" function
    remote.foo_user(function(err, val) {
      if(err) return t.fail("remote.foo_user call failed: " + err);

      t.equal(val, "bar_user");

      t.end();
    });
  });
});
