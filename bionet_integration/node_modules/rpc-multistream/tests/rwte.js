var fs = require('fs')
var rpc = require('../')
var test = require('tape')
var through = require('through2');
var tmp = require('tmp')

// tape test for rpc-multistream functionality:
// multiple streams:
//    read a file
//    write to a file
//    through uppercase function
//    error thrown

test('rwte', function (t) {
    var server = rpc({
        foo: rpc.syncReadStream(function() {
            return fs.createReadStream('tests/foo.txt', {encoding: 'utf8'});
        }),
        bar: rpc.syncWriteStream(function(filepath) {
            return fs.createWriteStream(filepath, {encoding: 'utf8'});
        }),
        baz: rpc.syncStream(function() {
            // create duplex transform stream with utf8 input and output
            return through({encoding: 'utf8', decodeStrings: false}, function(data) {
                this.push(data.toUpperCase());
            });
        }),
        bad: rpc.syncStream(function() {
            throw new Error("Something bad happened");
        })
    });
    t.plan(8)
    var client = rpc();
    client.pipe(server).pipe(client);
    client.on('methods', function(methods) {
        t.equal(typeof methods.foo, 'function', 'methods.foo is a function')
        t.equal(typeof methods.bar, 'function', 'methods.bar is a function')
        t.equal(typeof methods.baz, 'function', 'methods.baz is a function')
        t.equal(typeof methods.bad, 'function', 'methods.bad is a function')

        methods.foo().on('data', function(data) {
            t.equal(data,"I am the contents of foo.txt :)\n","foo creatReadStream foo.txt")
        })

        var tmpfile = tmp.fileSync()
        var outStream = methods.bar(tmpfile.name);
        outStream.write("Love!\n");
        outStream.end();
        fs.readFile(tmpfile.name, 'utf8', function(err,data) {
            if (err) t.fail(err)
            t.equal(data,"Love!\n","bar creatwritestream tmpfile")
        })

        var duplexStream = methods.baz()
        duplexStream.on('data', function(data) {
            t.equal(data,"EXONERATION\n","baz syncstream through")
        })
        duplexStream.write("exoneration\n");

        var bad = methods.bad()
        bad.on('error', function(err) {
            if (err.message == "Channel destroyed") return
            t.equal(err.message,"Something bad happened","bad syncstream throws error")
        })
    })
})
