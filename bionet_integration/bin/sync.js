#!/usr/bin/env node

var fs = require('fs');
var path = require('path');
var glob = require('glob');
var async = require('async');
var xtend = require('xtend');
var bc = require('bionet-client');

var settings = require(path.join(__dirname, '..', 'settings.js'));

var dataFilePath = path.join(__dirname, '..', '..', 'data', '*', '*.json');

function fail(err) {
  console.error(err);
  process.exit(1);
}

function saveVirtual(remote, data, cb) {
  var o = {
    freegenes: xtend(data, {}) // clone
  };

  if(data.virtual_id) {
    o.id = data.virtual_id;
    delete o.freegenes.virtual_id;
  }

  o.name = data.gene_name;

  o.description = "Submitted to the FreeGenes project";
  if(data.author && data.author.name) {
    o.description += " by " + data.author.name;
    if(data.author.affiliation) {
      o.description += " from " + data.author.affiliation;
    }
  }
  if(data.project_description) {
    o.description += " as part of the " + data.project_description + " project";
  } else if(data.author && data.author.project) {
    o.description += " as part of the " + data.author.project + " project";    
  }
  
  if(data.optimized_sequence) {
    o.sequence = data.optimized_sequence;
  }

  remote.saveVirtual(o, function(err, virtual) {
    if(err) return cb(err);

    if(!data.virtual_id) {
      data.virtual_id = virtual.id;
      cb(null, data);
    } else {
      cb();
    }

  });
}


function syncFile(remote, filepath, data, cb) {
  saveVirtual(remote, data, function(err, newData) {
    if(err) return cb(err);

    // the data didn't change so we don't need to write to file
    if(!newData) {
      return cb();
    }

    fs.writeFile(filepath, JSON.stringify(data, null, 2), {encoding: 'utf8'}, function(err) {
      if(err) return cb(err);

      cb(null, true);
    });
  });
}

function syncFiles(remote, cb) {

  var created = 0;
  var updated = 0;

  glob(dataFilePath, function(err, files) {
    if(err) return cb(err);

    if(!files || !files.length) {
      return cb(new Error("No files to synchronize. Expecting files here: " + dataFilePath));
    }

    console.log("Synchronizing " + files.length + " files...");
    
    async.eachSeries(files, function(file, cb) {
      
      fs.readFile(file, {encoding: 'utf8'}, function(err, data) {
        if(err) return cb(err);

        console.log("Synchronizing file");
        syncFile(remote, file, JSON.parse(data), function(err, didCreate) {
          if(err) return cb(err);

          if(didCreate) {
            created++;
          } else {
            updated++;
          }
          cb();
        });
      });

    }, function(err) {
      if(err) return cb(err);

      cb(null, created, updated);
    })
  });

}


bc(settings.node_url, {
  username: settings.username,
  password: settings.password
} ,function(err, done, remote) {
  if(err) return t.fail("failed to connect: " + err);

  syncFiles(remote, function(err, created, updated) {
    if(err) fail(err);


    console.log("Successfully created " + created + " and updated " + updated + " virtuals on " + settings.bionet_url)

    process.exit(0);
  });
});
