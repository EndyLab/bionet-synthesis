var fs = require('fs')
var rpc = require('../')
var test = require('tape')
var tmp = require('tmp')

// remote functions that return multiple streams as an array

test('syncArray', function (t) {
    tmp.setGracefulCleanup()
    var rwfile = tmp.fileSync()
    var server = rpc({
        multiRead: rpc.syncReadStream(function() { return [
            fs.createReadStream('tests/foo.txt', {encoding: 'utf8'}),
            fs.createReadStream('tests/foo2.txt', {encoding: 'utf8'})
        ]}, 2),
        multiRW: rpc.syncStream(function() { return [
            fs.createReadStream('tests/foo.txt', {encoding: 'utf8'}),
            fs.createWriteStream(rwfile.name, {encoding: 'utf8'})
        ]}, [{
            type: 'read',
            encoding: 'utf8',
            objectMode: false
        },{
            type: 'write',
            encoding: 'utf8',
            objectMode: false
        }])
    })
    var client = rpc()
    client.pipe(server).pipe(client)
    t.plan(4)
    client.on('methods', function(methods) {
        var streams = methods.multiRead()
        streams[0].on('data',function(data,err) {
            if (err) t.fail("mysterious error A: " + err)
            t.equal(data,"I am the contents of foo.txt :)\n","multiread foo")
        })
        streams[1].on('data',function(data,err) {
            if (err) t.fail("mysterious error B: " + err)
            t.equal(data,"I am the contents of foo2.txt :)\n","multiread foo2")
        })

        var rwStreams = methods.multiRW()
        rwStreams[0].on('data',function(data,err) {
            if (err) t.fail("mysterious error C: " + err)
            t.equal(data,"I am the contents of foo.txt :)\n","multirw foo")
        })
        rwStreams[1].write("I am the contents of rwfile")
        fs.readFile(rwfile.name,'utf8',function(err,data) {
            if (err) t.fail("mysterious error D: " + err)
            t.equal(data,"I am the contents of rwfile","multirw tmpfile")
        })
    })
})
