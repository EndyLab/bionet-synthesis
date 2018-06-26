
var bc = require('../index.js');
var test = require('tape')


test('simple', function(t) {

  t.plan(1);

  var url = "http://localhost:8000";

  bc(url, function(err, done, remote) {
    if(err) return t.fail("failed to connect: " + err);


    remote.foo(function(err, val) {
      if(err) return t.fail("remote.foo call failed: " + err);

      t.equal(val, "bar");

      t.end();
      process.exit(0);
    });
  });
});
