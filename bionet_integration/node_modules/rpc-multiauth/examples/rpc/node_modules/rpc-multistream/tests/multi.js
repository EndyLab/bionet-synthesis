var fs = require('fs')
var rpc = require('../')
var test = require('tape')

// tape test for rpc-multistream functionality:
// calling multiple functions asynchronously

test('multi', function (t) {

    t.plan(6)

    var server = rpc({
        foo: function(str, cb) {
            str = str.toUpperCase();
            cb(null, str);
        },
        bar: function(path, cb) {
            var stream = fs.createReadStream(path);
            cb(null, stream);
        },
        baz: rpc.syncReadStream(function(path) {
            return fs.createReadStream(path);
        })
    });

    var client = rpc();

    client.pipe(server).pipe(client);

    client.on('methods', function(methods) {

        t.equal(typeof methods.foo, 'function', 'methods.foo is a function')
        t.equal(typeof methods.bar, 'function', 'methods.bar is a function')
        t.equal(typeof methods.baz, 'function', 'methods.baz is a function')

        methods.foo("rose", function(err, msg) {
            t.equal(msg, "ROSE", "foo toUpperCase ROSE")
        });

        methods.bar("tests/foo.txt", function(err, stream) {
            stream.setEncoding('utf8')
            stream.on('data', function(data) {
                t.equal(data,"I am the contents of foo.txt :)\n","bar creatReadStream foo.txt")
            })
        });

        methods.baz("tests/foo.txt").on('data', function(data) {
            t.equal(data,"I am the contents of foo.txt :)\n","baz creatReadStream foo.txt")
        })

    });

})
