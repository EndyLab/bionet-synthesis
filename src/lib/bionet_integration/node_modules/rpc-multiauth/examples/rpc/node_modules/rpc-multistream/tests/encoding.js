var from = require('from2')
var fs = require('fs')
var rpc = require('../')
var test = require('tape')

// using both binary, string and object streams

test('encoding', function (t) {
    var server = rpc({
        foo: rpc.syncReadStream(function() {
            return fs.createReadStream('tests/foo.txt')
        }, {
            encoding: null, // binary
            objectMode: false
        }),
        bar: rpc.syncReadStream(function() {
            return fs.createReadStream('tests/foo.txt')
        }, {
            encoding: 'utf8',
            objectMode: false
        }),
        baz: rpc.syncReadStream(function() {
            var i = 0
            return from.obj(function(size, next) {
                if (i++) return next(null, null)
                next(null, {
                    hoopy: 'frood'
                })
            })
        }, {
            objectMode: true
        })
    }, {
        objectMode: false
    })
    var client = rpc(undefined, {
        objectMode: false,
        explicit: true
    })
    client.pipe(server).pipe(client)
    t.plan(3)
    client.on('methods', function(methods) {
        methods.foo().on('data', function(data) {
            t.deepEqual(data,new Buffer("I am the contents of foo.txt :)\n", 'binary'),"foo binaryStream")
        })
        methods.bar().on('data', function(data) {
            t.equal(data,"I am the contents of foo.txt :)\n","bar utf8Stream")
        })
        methods.baz().on('data', function(data) {
            t.deepEqual(data,{ hoopy: 'frood' },"baz objStream")
        })
    })
})
