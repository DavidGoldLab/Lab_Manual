<!-- #region -->
# Lesson 3: Git and GitHub

© David Gold. Except where the source is noted, this work is licensed under a [Creative Commons Attribution License CC-BY 4.0](https://creativecommons.org/licenses/by/4.0/).

## 3.1. Installing Git

The Gold lab uses Git for version control when sharing files, data, and code.

Put briefly, Git is a __version control system__ that allows multiple people to collaborate on the same code while keeping track of the changes each person makes. Git allows locally stored __repositories__ of code to be synchronized against a remotely stored master copy.

Install the latest version of Git with Homebrew:

```
brew install git
```

## 3.2. Signing up for a GitHub account

GitHub is a company (owned by Microsoft) that provides hosting for software development. Github connects with Git so that local computers can easily interact with code stored on the website.

Go to [http://github.com/](http://github.com/) to get an account. You should register with your UC Davis email address (you get free repositories as an academic).

Yellow boxes (alert-warning)
<div class="alert alert-block alert-warning">
<b>Warning: Be thoughtful about the username you choose.</b> This website is part of your professional profile--for some careers your GitHub page will be more important than your resume.
</div>


## 3.3. Setting your username in Git

You can change the name that is associated with your `git` commands using the `git config` command. The name you set will be visible in any future changes (also called "__commits__") you push to GitHub from Terminal. You do not have to use your real name if you don't want to.

Configure your Github username. __Don't forget to replace "Barbara McClintock" with your prefered name!__:

```
git config --global user.name "Barbara McClintock"
```

Now configure your Github email. __Don't forget to replace "your_email@example.com" with your actual email address!__

```
git config --global user.email "your_email@example.com"
```

Finally, set up BBEdit as your prefered text editor:

```
git config --global core.editor "bbedit -w"
```

## 3.4. Forking a repository

GitHub is full of repositories with useful pieces of code. You might not be a collaborator on a repository, but you want to play with the code and perhaps make modificiations to it. You also want to keep the code up to date in case the author(s) make any changes.

An easy way to do this is through __forking__. To fork a repository on GitHub, go to the repository on Github and click the "fork" button. __Be sure you are logged into GitHub when you do this__. 

Click [this link](https://github.com/DavidGoldLab/Gold_Lab_Training) to go to the GitHub page for our lab.

The fork button is in the upper right. Click the button and you  have a __fork__ of our repository on your GitHub account!

## 3.5. Cloning your fork to your local machine

Now that you have a forked copy of the repository associated with your GitHub account, you can clone it onto your computer. To be clear, the __local version__ of your repository is stored on your computer's hard drive; it is distinct from the repository stored on GitHub (although the two can interact, as you will see soon).

Let's store all of your repositories in a directory called `git` in your home directory. We can make a directory by typing the following command into Terminal:

```
mkdir ~/git
```

Now you can clone your forked copy of the repository onto your local machine. To do this, navigate your browser to the forked copy of the `Gold_Lab_Training` repository on your account (this is where clicking the "Fork" button took you in your browser). The browser URL will be: `https://github.com/YOUR_USERNAME/Gold_Lab_Training`, and the top left of the website will say something like "`yourusername/Gold_Lab_Training` forked from `DavidGoldLab/Gold_Lab_Training`".


Now you can clone it with the following command in Terminal. __Make sure to replace `https://github.com/YOUR_USERNAME/Gold_Lab_Training` in the code below with your personal URL:

```
cd ~/git
git clone https://github.com/YOUR_USERNAME/Gold_Lab_Training
```

You now have a local copy of your own fork of the bootcamp repository. You can add files and edit it. When you commit and push, it will all be on your account, and the master repository will not see the changes.

## 3.6. Syncing your forked repository to the upstream repository

You want to be able to sync your forked repository with the original repository (also called the __upstream repository__) in case there are any updates. What we now need to do is link your local copy of the repository to the upstream Gold_Lab_Training repository.

Using Terminal, navigate to your local copy of the Gold_Lab_Training folder. You can then use the `git remote` command to see which remote (online) repositories are associated with this repository. Adding the `-v` (__verbose__) flag will provide additional information, including the URLs associated with the remote repositories:

```
cd Gold_Lab_Training
git remote -v
```

Part of your output should look something like this:

    origin	https://github.com/YOUR_USERNAME/Gold_Lab_Training (fetch)
    origin	https://github.com/YOUR_USERNAME/Gold_Lab_Training (push)

This shows you the URL to your repository. Now you can connect your version of the repository to the upstream version:

```
git remote add upstream https://github.com/DavidGoldLab/Gold_Lab_Training
```

Now try doing `git remote -v`, and you will see that you are now connected to the upstream repository.

    origin	https://github.com/YOUR_USERNAME/Gold_Lab_Training (fetch)
    origin	https://github.com/YOUR_USERNAME/Gold_Lab_Training (push)
    upstream	https://github.com/DavidGoldLab/Gold_Lab_Training (fetch)
    upstream	https://github.com/DavidGoldLab/Gold_Lab_Training (push)


```{tip}
If you ever want to delete a repository, go to the git directory and remove the repository using the force (`-f`), and recursive (`-r`) flags:

cd ~\git\
rm -rf Gold_Lab_Training
```
    
## 3.7. Making changes in the local repository folder

If you want to change file names for annything in your foler, the  `mv` command with `git`. That way, Git will keep track of the naming changes you made:

    git mv [Original_File/Path] [New_File/Path] 

## 3.8. Communicating between the local repository and GitHub

Now that you've got a local copy and a copy on your GitHub account, there are five commands you need to know in order to interact with your forked repository on Github:

- __Add__ - If you add new files to your local repository, you need to use this command before using `commit` (see below).
- __Commit__ - This command records changes you have made to the repository. Think of it as a snapshot of the current status of the project. Commits are done locally (in other words, on your personal computer).
- __Push__ - This command "pushes" the recent commit from your local repository up to GitHub. If you're the only one working on a repository, pushing is fairly simple. If there are others accessing the repository, you may need to pull before you can push.
- __Pull__ - This command "pulls" any changes from the GitHub repository and merges them into your local repository.

## 3.9. Saving changes in local Git repository to Github

Make sure you are in the relevant repository directory in terminal. Then `add` the changes you made using git:

```
git add --all
```

`Commit` the changes. __Replace 'commit_message' with a short message detailing the changes you made (this is useful for version control)__:

```
git commit -m 'commit_message'
```

Copy files from your local repository to GitHub using `push`

```
git push
```

You may be asked to provide your GitHub username and password

<!-- #endregion -->
