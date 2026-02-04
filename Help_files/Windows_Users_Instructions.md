Windows users can run the LSPCA on Docker Desktop by the following steps.

1. Install [Docker Desktop](https://www.docker.com/products/docker-desktop/) for Windows and create an account. 

2. Next you would need to downlowad the LSPCA package from the Guithub repository "jamnamdari\LSPCA" and extract the zipped file.

3. Open Windows PowerShell and cd to the LSPCA folder.

4. run the following commands in Windows PowerShell

```
docker build . -t lspca
```
and 

```
 docker run --rm -d -p 8787:8787 -e USER=test -e PASSWORD=1234 -e ROOT=TRUE -v "${PWD}:/home/test/project" --name lspca lspca
```

5. type "http://localhost:8787" in chrome browser
6. login with
   
   username: test
   
   password: 1234

7. loading all package dependencies by running: 
```
source("project/LSPCA/R/setup.R")
```

LSPCA package is loaded and ready to use!


