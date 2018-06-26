var fs = require('fs')
var path = require('path')
var rpc = require('../')
var test = require('tape')
var tmp = require('tmp')
var through = require('through2')

// asynchronous remote functions that return streams

test('asyncStream', function (t) {
    tmp.setGracefulCleanup()
    var foofile = tmp.fileSync()
    var bazfile = tmp.fileSync()
    var server = rpc({
        foo: function(cb) {
            var w = fs.createWriteStream(foofile.name, {encoding: 'utf8'})
            cb(null, w)
        },
        bar: function(cb) {
            var r = fs.createReadStream('tests/foo.txt', {encoding: 'utf8'})
            cb(null, "bar says hi", r)
        },
        baz: function(cb) {
            var w = fs.createWriteStream(bazfile.name, {encoding: 'utf8'})
            var r = fs.createReadStream('tests/foo.txt', {encoding: 'utf8'})
            cb(null, r, w, {chiao: "hiya!"})
        },
        duper: function(cb) {
            var ds = through({encoding: 'utf8', decodeStrings: false}, function(data) {
                data = data.toUpperCase()
                this.push(data)
            })
            cb(null, ds)
        }
    })
    var client = rpc()
    client.pipe(server).pipe(client)
    t.plan(9)
    client.on('methods', function(methods) {
        methods.foo(function(err, w) {
            t.pass("foo")
            w.write("woop!")
            w.end()
            fs.readFile(foofile.name,'utf8',function(err,data) {
                if (err) t.fail("mysterious error A: " + err)
                t.equal(data,"woop!","foo file")
            })
        })
        methods.bar(function(err, msg, r) {
            t.equal(msg,"bar says hi","bar")
            r.on('data', function(data) {
                t.equal(data,"I am the contents of foo.txt :)\n","bar foo.txt")
            })
        })
        methods.baz(function(err, r, w, msg) {
            t.deepEqual(msg,{chiao:'hiya!'},"baz chiao")
            r.on('data', function(data) {
                t.equal(data,"I am the contents of foo.txt :)\n","baz foo.txt")
            })
            w.write("a treat for your tummy!")
            w.end()
            fs.readFile(bazfile.name,'utf8',function(err,data) {
                if (err) t.fail("mysterious error B: " + err)
                t.equal(data,"a treat for your tummy!","baz file")
            })
        })
        methods.duper(function(err, dupStream) {
            t.pass("duper returned")
            dupStream.on('data',function(data) {
                t.equal(data,"MAKING THIS UPPERCAAAAASE","duper uppercase")
            })
            dupStream.write("making this uppercaaaaase")
        })
    })
})
