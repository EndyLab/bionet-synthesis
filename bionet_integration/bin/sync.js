#!/usr/bin/env node

var fs = require('fs');
var path = require('path');
var glob = require('glob');
var async = require('async');
var xtend = require('xtend');
var bc = require('bionet-client');

// The maximum file version supported by this program.
// Any file with no version field is assumed to be version '0.0.0'.
// Any file with a version above this will be skipped
MAX_SUPPORTED_VERSION = '0.0.0';

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
  if(data.author && data.author.project) {
    o.description += " as part of the project: " + data.author.project;
  } else if(data.project_description) {
    o.description += " as part of the project: " + data.project_description;
  } else if(data.description) {
    o.description += " as part of the project: " + data.project_description;
  }

  if(data.optimized_sequence) {
    o.sequence = data.optimized_sequence;
  } else if(data.sequence && data.sequence.optimized_sequence) {
    o.sequence = data.sequence.optimized_sequence;
  }

  remote.saveVirtual(o, function(err, id) {
    if(err) return cb(err);

    if(!data.virtual_id) {
      data.virtual_id = id;
      cb(null, data);
    } else {
      cb();
    }

  });
}

// parses a version like 1.12.3 into a number like 10000012000003
function parseVersion(str) {
  if(!str) throw new Error("Failed to parse version string: " + str);

  var fields = str.split('.');
  if(fields.length != 3) {
    throw new Error("Failed to parse version string: " + str);
  }

  var n;
  fields = fields.map(function(s) {
    n = parseInt(s);
    if(n > 999999) {
      throw new Error("Failed to parse version string: " + str);
    }
    return n;
  });

  return fields[0] * 1000000000000 + fields[1] * 1000000 + fields[2];
}

function syncFile(remote, filepath, data, cb) {
  if(data.version) {
    try {
      if(parseVersion(data.version) > parseVersion(MAX_SUPPORTED_VERSION)) {
        console.warn("Encountered unsupported file version:", data.version);
        return cb(null, false, true); // skip
      }
    } catch(err) {
      return cb(err);
    }
  }

  if(data.status && !data.status.will_build) {
    return cb(null, false, true); // skip
  }

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
  var skipped = 0;

  glob(dataFilePath, function(err, files) {
    if(err) return cb(err);

    if(!files || !files.length) {
      return cb(new Error("No files to synchronize. Expecting files here: " + dataFilePath));
    }

    console.log("Synchronizing " + files.length + " files...");
    
    async.eachSeries(files, function(file, cb) {
      
      fs.readFile(file, {encoding: 'utf8'}, function(err, data) {
        if(err) return cb(err);

        console.log("Synchronizing file:", file);

        syncFile(remote, file, JSON.parse(data), function(err, didCreate, skipped) {
          if(err) return cb(err);

          if(skipped) {
            skipped++;
          } else {
            if(didCreate) {
              created++;
            } else {
              updated++;
            }
          }
          cb();
        });
      });

    }, function(err) {
      if(err) return cb(err);

      cb(null, created, updated, skipped);
    })
  });

}


bc(settings.node_url, {
  username: settings.username,
  password: settings.password
} ,function(err, done, remote) {
  if(err) fail("failed to connect: " + err);

  syncFiles(remote, function(err, created, updated, skipped) {
    if(err) fail(err);


    console.log("Successfully skipped " + skipped + ", created " + created + " and updated " + updated + " virtuals on " + settings.node_url)

    process.exit(0);
  });
});
