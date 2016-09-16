#!/bin/bash

# enter subtree names with space between each entry
subtrees=(python R)

echo "Number of subtrees: ${#subtrees[*]}"

# Push each subtree with standard command
echo "Pushing with 'git subtree push' to master branch using --squash"
for item in ${subtrees[*]}
do
    git subtree push --prefix=$item --squash $item master
    echo ""
done
