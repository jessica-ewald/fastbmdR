## Testing
If changes are made to the fastbmdR_main.R or fastbmdR_utils.R, run `Rscript test/test.R` from the root fastbmdR directory.


## Documentation
As long as Roxygen comments and tags are included above all functions that should be callable (exported), we can automatically generate the necessary NAMESPACE and man/*.md pages with:
```r
setwd('path-to-fastbmdR-repo')
devtools::document()
```
This ensures that functions can be called and that documentation will show up when someone enters ?function-name in R. 


## Releases
After commiting changes, enter:
```
git tag v#.#.#
git push origin v#.#.#
```

Then, go to the Github interface, go to Releases, and find the tag that you just pushed. The release should be in draft form. Finish adding whatever notes are needed and publish the release. In the case that you immediately notice a bug, we can delete the release and republish. This is bad practice, but since fastbmdR has approximately 3 users, it should be fine.

```
git tag -d v1.0.0
git push origin :refs/tags/v1.0.0  # Delete remote tag

# Make changes ...
git add .
git commit -m "Fix bug in v1.0.0"

# Recreate the tag
git tag v1.0.0
git push origin v1.0.0
```
