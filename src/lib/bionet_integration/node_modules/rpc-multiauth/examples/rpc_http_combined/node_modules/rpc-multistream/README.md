
rpc-multistream is similar to [rpc-stream](https://github.com/dominictarr/rpc-stream) but has support for remote functions that return streams like [multilevel](https://github.com/juliangruber/multilevel) and works without telling the client which remote functions are available like [dnode](https://github.com/substack/dnode). 

rpc-multistream uses [mux-demux](https://github.com/dominictarr/mux-demux), a slighty tweaked version of [rpc-stream](https://github.com/juul/rpc-stream) and borrows some code from [multilevel](https://github.com/juliangruber/multilevel).


# Usage #

```
var fs = require('fs');
var rpc = require('rpc-multistream');

var server = rpc({
    foo: rpc.readable(function() {
        return fs.createReadStream('foo.txt');
    }),
    bar: function(cb) {
        console.log("bar called");
        cb(null, "bar says hi");
    }
});

var client = rpc();

client.pipe(server).pipe(client)

client.on('remote', function(remote) {

    var stream = remote.foo();
    stream.on('data', function(data) {
        console.log(data);
    });

    remote.bar(function(err, msg) {
        console.log(msg);
    });
});
```

# License #

MIT

Copyright (c) 2014 Marc Juul <juul@sudomesh.org>
Copyright (c) 2013 Julian Gruber <julian@juliangruber.com>

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
