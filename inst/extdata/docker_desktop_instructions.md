## Docker Desktop Instructions 

1. Install Docker Desktop application from [https://www.docker.com/products/docker-desktop/](https://www.docker.com/products/docker-desktop/).

2. Pull the docker image from Docker Hub using the following command on the terminal:

```
docker pull amanukyan1385/rstudio-voltron:main
```

3. Open Docker Deskop app and go to `Images`

4. Find `amanukyan1385/rstudio-voltron:main` and click `Run` on the right. 

5. Enter `8787` to Ports, `PASSWORD` to `Variable` and `yourpassword` (you can also give a custom password) to `Value` under Environment variables. If you want to add folders to the docker image, you can update `Volumes` by entering the path to the folder in your local machine in `Host path`, and then by entering the path in the docker container under `Container path` (e.g. /home/rstudio/workshop/).  

6. Click `Run`

7. Start the RStudio session from the browser at `http://localhost:8787/` and enter `rstudio` as username and `yourpassword` as password. 
