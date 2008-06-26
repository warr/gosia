#!/bin/bash

list=$(git tag -l |grep ^200)

for t in $list ; do
    echo "Creating archive for gosia $t"
    git archive --format=tar --prefix=gosia_$t/ $t | gzip -9 > /tmp/gosia_$t.tgz
done
