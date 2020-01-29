---
title: "Ecological Genomics Notebook"
author: "Anoob Prakash"
date: "22/01/2020"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

### **Affiliation**:  Keller Lab; Department of Plant Biology, University of Vermont  
### **E-mail contact**:  anoob.prakash@uvm.edu  
### **Start Date**: 2020-01-13  
### **End Date**: 2020-05-08  
### **Project Descriptions**: Climate change adaptation in red spruce  

**************************************************************************

### Table of Contents:  

- [ ] [Working with fastq](./fastq) 

___________________________________________________________________________  

## General notes  
**Log into UNIX server through RStudio terminal**  
`require(ssh)`
`session <- ssh_connect("ap1@pbio381.uvm.edu")`
`print(session)`

Go to terminal on RStudio, type in the following to get root access:  

`/usr/bin/env bash -l`  

type in your password  

You have root access now:
 type in the details to login to the server  
 
 `ssh servername`
 
 For this course, type in netID followed by the server ID. For example:
 
 `ssh ap1@pbio381.uvm.edu`  

 Next type in your net ID password, and now you can accesss the UNIX shell like you would be accesssing through the terminal on MAC or through PUTTY on Windows. Additonal tip for Windows users, there is an alterantive to putty on windows store. Just download and install ubuntu and it works as a linux shell. It seems more integrated and part of the windows system, however you cannot save login credentials like you can on putty.  
   
 
 More details on RStudio terminal can be  found here: https://support.rstudio.com/hc/en-us/articles/115010737148-Using-the-RStudio-Terminal   
 
***
**Tutorial on bash :** https://pespenilab.github.io/Ecological-Genomics/Tutorial/2020-01-22_Command_Line_Unix.html   

**Course website:** https://pespenilab.github.io/Ecological-Genomics/  

***  
