var rpc = require('../')
var test = require('tape')

// Simple rpc-multistream test of asynchronous calling

test('toUpperCase', function (t) {

    t.plan(1)

    var server = rpc({
        foo: function(str, cb) {
            str = str.toUpperCase();
            cb(null, str);
        }
    });

    var client = rpc();

    client.pipe(server).pipe(client);

    client.on('methods', function(methods) {
        methods.foo("rose", function(err, msg) {
            t.equal(msg, "ROSE", "toUpperCase")
        });
    });

})
