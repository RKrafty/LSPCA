Mac users would need to install [Docker Desktop](https://www.docker.com/products/docker-desktop/) to be able to run the R package. Instructions on how to install Docker and run the LSPCA package are provided below. Note that Docker Desktop is supported on the current and two previous major macOS releases.

1. Install brew (if not already installed):
```
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
```
2. Setup Docker Desktop:
```
brew install --cask docker
```

Next you wowuld need to downlowad the LSPCA package from the Guithub repository "jamnamdari\LSPCA" 

3. cd into LSPCA project and run:
```
make
```

4. type "http://localhost:8787" in chrome browser
5. login with
   
   username: test
   
   password: 1234

7. loading all package dependencies by running: 
```
source("project/R/setup.R")
```

LSPCA package is loaded and ready to use!


