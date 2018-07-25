#!/bin/sh
mess="Update documentation"
if [ $1 ] ; then
    mess=$1
fi
git checkout master
nimble gendoc
git checkout gh-pages
cp -r docs/html/* .
rm -r *.idx
git add *
git commit -m "$mess"
