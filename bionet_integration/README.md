
This is a small node.js command-line program to synchronize data between the freegenes project and a bionet node.

# Dependencies

Ensure that you have a recent version of node.js installed. Then in this directory run:

```
npm install
```

# Setup

Copy `settings.js.example` to `settings.js`. Ensure that the `settings.js` file is only readable by those who should have access to the bionet node account credentials, e.g:

```
touch settings.js
chmod 600 settings.js
cat settings.js.example > settings.js
```

Now edit `settings.js`. You will need to set the bionet node url and user credientials. 

# Usage

When running `./bin/sync.js` the `.json` files in the `../data` directory are read and created as virtuals on the bionet node, or updated if they already exist on the bionet node. If a new virtual is created for a `.json` file then the `.json` file will be updated to include the `virtual_id` used to identify the associated virtual on the bionet node.

A normal workflow will be something like:

```
git pull
./bin/sync.js
git commit -a -m "synchronized with bionet node"
git push
```

# Testing

To test that the connection to the server is working, without modifying any data, run:

```
npm run test
```

# Copyright and license

Copyright 2018 BioBricks Foundation

License: AGPLv3
