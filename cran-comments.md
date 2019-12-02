## Test environments
* local R installation, R 3.6.1
* ubuntu 16.04 (on travis-ci), R 3.6.1
* win-builder (devel)

## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new release.




## Additional CRAN comments (Dec 1, 2019)

### Please make sure that you do not change the user's options, par or working directory. If you really have to do so, please ensure with an *immediate* call of on.exit() that the settings are reset when the function is exited.

I have used the `on.exit()` function in the `plotInter()`, `plot.sail()` function to fix this issue

#### Do really all of your examples need to be wrapped in if(interactive()) ?

I have removed the `interactive()` wrapper around most of the examples. 


## Additional CRAN comments (Nov 25, 2019)

### Please only capitalize names and sentence beginnings in the description.

Fixed

### Please ensure that you do not use more than 2 cores in your examples.

Fixed

### \dontrun{} should be only used if the example really cannot be executed (e.g. because of missing additional software, missing API keys, ...) by the user. That's why wrapping examples in \dontrun{} adds the comment ("# Not run:") as a warning for the user. Does not seem necessary, as you wrapped it in if(interactive())

I have removed \dontrun{} from all of the examples. 


### You write information messages to the console that cannot be easily suppressed. It is more R like to generate objects that can be used to extract the information a user is interested in, and then print() that object. Instead of print()/cat() rather use message()/warning()  or if(verbose)cat(..) if you really have to write text to the console. (except for print() and summary() functions)

I have replaced `cat()` in the `kkt()` function with `warning()`


#### Please halve these tests runtimes: The overall threshold for the whole package check on CRAN is roughly 10 min. This had 15 min.

Fixed. I have significantly reduced the tests times. 

