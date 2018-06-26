var from = require('from2');
var fs = require('fs')
var rpc = require('../')
var test = require('tape-catch')

// tape test for rpc-multistream functionality:
// opts.explicit

test('opts.explicit.on', function (t) {
    var server = rpc({
        foo: rpc.syncReadStream(function() {
            return fs.createReadStream('tests/foo.txt', {encoding: 'utf8'});
        }),
        bar: rpc.syncReadStream(function() {
            var i = 0;
            return from.obj(function(size, next) {
                if(i++) return next(null, null);
                next(null, {
                    hoopy: 'frood'
                });
            });
        }, { objectMode: true }) // bar is objectMode
    }, { objectMode: false }); // foo is not
    t.plan(5)
    var client = rpc(undefined, {
        objectMode: true, // different default than the server
        explicit: true // with explicit set it doesn't matter if defaults differ
    });
    client.pipe(server).pipe(client);
    client.on('methods', function(methods) {
        t.equal(typeof methods.foo, 'function', 'explicit client: methods.foo is a function')
        t.equal(typeof methods.bar, 'function', 'explicit client: methods.bar is a function')

        methods.foo().on('data', function(data) {
            t.equal(data,"I am the contents of foo.txt :)\n","explicit client: foo creatReadStream foo.txt")
        })

        var objStream = methods.bar();
        objStream.on('data', function(data) {
            t.equal(typeof data,"object","explicit client: bar gives an object")
            t.deepEqual(data,{hoopy: 'frood'},"explicit client: bar object is correct")
        });
    })
})
test('opts.explicit.off', function (t) {
    var server = rpc({
        foo: rpc.syncReadStream(function() {
            return fs.createReadStream('tests/foo.txt', {encoding: 'utf8'});
        }),
        bar: rpc.syncReadStream(function() {
            var i = 0;
            return from.obj(function(size, next) {
                if(i++) return next(null, null);
                next(null, {
                    hoopy: 'frood'
                });
            });
        }, { objectMode: true }) // bar is objectMode
    }, { objectMode: false }); // foo is not
    t.plan(3)
    var client = rpc(undefined, {
        objectMode: true, // different default than the server
        explicit: false // this causes errors, now we have to catch them
    });
    client.pipe(server).pipe(client);
    client.on('methods', function(methods) {
        t.equal(typeof methods.foo, 'function', 'explicit client: methods.foo is a function')
        t.equal(typeof methods.bar, 'function', 'explicit client: methods.bar is a function')

        methods.foo().on('error', function(err) {
            t.equal(err.message,"Could not parse row I am the contents of foo.txt :)...","stringstream error")
        })

        var objStream = methods.bar();
        objStream.on('data', function(data) {
            t.equal(typeof data,"object","explicit client: bar gives an object")
            t.deepEqual(data,{hoopy: 'frood'},"explicit client: bar object is correct")
        });
    })
})
