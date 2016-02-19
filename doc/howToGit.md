First time setup:
  About you:
  ```
  git config --global user.name yourUserName
  git config --global user.email yourEmail@foo.bar
  ```
  Set up your favourite editor to be used:
  ```
  git config --global core.editor <editorname>
  ```
  with editorname being whatever you prefer (vim, emacs, nano...)

  
  Setting up a diff tool:
  ```
  git config --global diff.tool <toolname>
  ```
  Again with your choice of the difftool (tkdiff, meld...)

General use:

  get the repository the first time:
  ---------------------------------
  ```
  git clone https://github.com/abieler/dsmc.git
  ```


  modifying the file fileName.jl and commiting changes:
  ---
  ```
  git add fileName.jl
  git commit
  git push
  ```
  
  
* if you are not happy with the changes you made to the file, you can again set it back to 
  whatever the status is on the master by checking out the corresponding file from there with:
  ```
  git checkout filename
  ```


* creating a new branch:
  ```
  git checkout -b newBranchName
  ```
  (then do the git add and commits to work in 
  the new branch as usual)


* push your new branch to remote:
  ```
  git push --set-upstream origin newBranchName
  ```

* get branch from remote to your local box
  ```
  git fetch origin branchName
  git checkout branchName
  ```
  

* changing between branches:
  ```
  git checkout branchName
  git checkout master (to check back to master)
  ```


* merging new branch with master:
  ```
  git checkout master (change back to master branch)
  git merge newBranchName
  git branch -d newBranchName (delete branch)
  ```
