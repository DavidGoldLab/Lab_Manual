# Lesson 1: Basic Command Line Skills

© David Gold. Except where the source is noted, this work is licensed under a [Creative Commons Attribution License CC-BY 4.0](https://creativecommons.org/licenses/by/4.0/).

## 1.1. Getting started with the command line

Now, as you did in lesson 1, open up your Terminal application

## 1.2. Navigating around directories

When you click around the folders and subfolders of your computer, you are actually navigating __directories__. Directories are a hierarchical file system your computer uses to store and organize files.

Let's start by using the `pwd` __(print working directory)__ command to figure out what directory we're in. Type the following into terminal:

````
pwd
````

You should get something like:

    /Users/{Your_Computer_Username}

`pwd` tells you the __path__ of your current directory. The path shows where a file or folder is stored on your computer. The path lists all of the parent directories in the hierarchy, separated by slashes (`/`), all the way up to the root directory, which is signified by the initial `/`. When you start up Terminal you should be in your __home directory__.

To list all files and folders in the current directory, we employ the `ls` __(list)__ command:

````
ls
````

## 1.3. Changing directories

Let's move arond the computer using the `cd` __(change directory)__ command. If you are in your home directory, you should se the "Desktop" folder when you use the `ls` command. This is because the "Desktop" is a subfolder in your home directory (in other words, your home directory is the __parent directory__ and your desktop is the __child directory__).

You can navigate to the Desktop using the `cd` command:

````
cd Desktop
````

```{Tip}
Instead of typing out the whole name to a directory or file, try typing part of the name out and hitting the `tab` button. Assuming there are not multiple files and/or folders with the same start, Terminal will auto-fill the word based on the objects in your directory. 

In this case, typing "D" + "tab" will probably not work, because you probably have other folders that start with "D" (e.g. Downloads, Documents). But typing "D" + "e" + "tab" probably will work. You should get comfortable using tab to auto-fill names; it makes command line work much easier!
```

Use the `ls` command again and you can see all of the files and folders currently on your desktop.

Another way to specify your home directory is by its shortcut, `~`. So you can always go direclty to your home directory using the following command:

````
cd ~
````

If you know the path to a folder, you can go directly to it using the `cd` command:

````
cd ~/Desktop
````

One last trick; you can use periods as a shortcut to move around directories. One period (`.`) signifies the folder you are in. Two periods (`..`) signifies the parent folder. So if you wanted to navigate up two directories you could do the following:

````
cd ../../
````

## 1.4. View your working directory in Finder

Here is a nice trick if you want to "__open__" the directory in a more traditional "Finder" window:

````
open .
````

The period means that you want to open the current directory. You could open a different folder by providing its path.


## 1.5. Making and removing directories

Let's start this exercise by making sure we are in the home directory:

````
cd ~
````

### 1.5.1. Making directories

To make a directory, the command is `mkdir` (__make directory__) followed by the name of the directory you want to create. For example, to make a directory called `test` type the following into Terminal:

````
mkdir test
````

You now have an empty directory called test. You can see it if you list the contents of your current working directory.

````
ls
````

You can move into this directory with the `cd` command.

````
cd test
````

Verify that you are in the correct directory by checking your path:

````
pwd
````

In response the computer should report a path that ends in "test"; for example:

    /Users/davidgold/test

Let's now move back into the parent directory:

````
cd ../
````

In the same way that one dot (`.`) represents the current directory (as mentioned in the `open .` command), two dots (`..`) means one directory up.

### 1.5.2. Removing (empty) directories

We do not need (nor want) this test directory, so let's delete it. To delete an *empty* directory, the command is `rmdir` (__remove directory).

````
rmdir test
````

To delete directories with files in them, you need the `rm` command (discussed later).

## 1.6. Creating and viewing text files

For the next exercise, I want to start by making a folder and adding a text file to it:

Make the "test2" directory with `mkdir`:

````
mkdir test2
````

Move into the "test2" directory with `cd`:

````
cd test2
````

### 1.6.1. Making a text file with nano

Now we're going to use a new command called `nano`. __Nano__ is a simple text-editor that you can acess from Terminal. You can open the text-editor by simply typing `nano` into terminal, or you can provide a filename for the text document you want to create:

````
nano Textfile.txt
````

This will open the Nano text-editor in Terminal. You can add any text you want, here's an example:
    
    Hello world
    
Once you've written some text, use the `command` + `x` keys to exit Nano. Nano will ask the following:

`Save modified buffer (ANSWERING "No" WILL DESTROY CHANGES) ?`

Press the `y` key to save your file. Nano will then double-check what you want to name the file:

`File Name to Write: Textfile.txt`

If you never provided a file name, you will have the opporutnity to do so now. If you are happy with the file name, hit the `enter / return` key. This will take you back to the normal Terminal window.

Verfiy that you now have a text file called "Textfile.txt" in your Test folder:

````
ls
````

You should get "Textfile.txt" in response.

### 1.6.2. Viewing text files

You can use terminal to look at the contents of a text file without opening it. Copy the command below and paste it into Terminal to create a new file called "File2.txt"

```
echo -e '1\n2\n3\n4\n5\n6\n7\n8\n9\n10\n11\n12\n13\n14\n15\n16\n17\n18\n19\n20' > File2.txt
```
File2.txt is 20 lines long, containing one number on each line.

To quickly see the first ten lines of this file, you can use the `head` command:

```
head File2.txt
```

It will report the first ten lines, which contain the numbers 1 through 10:

    1
    2
    3
    4
    5
    6
    7
    8
    9
    10

To quickly see the last ten lines of this file, you can use the `tail` command:

```
tail File2.txt
```

It will report the last ten lines, which contain the numbers 11 through 20:

    11
    12
    13
    14
    15
    16
    17
    18
    19
    20

By defauly, `head` and `tail` report ten lines. You can change that number with the `-n` flag:

```
head -n 18 File2.txt
```

The above command will report the first 18 lines of the document.

## 1.7. Moving and renaming files

Let's return to our original text file ("Textfile.txt"). You can rename it using the `mv` (__move__) command.

````
mv Textfile.txt Newfile.txt
````

Use the `ls` command again and you should see "Newfile.txt" instead of "Textfile.txt"

You can also move the file (with or without changing the file name) using this command.:

````
mv Newfile.txt ~/Desktop/Newfile.txt
````

## 1.8. Removing files and directories

You have learned some good skills, but now you have some junk spread all over your computer, including a text file on your desktop ("~/Desktop/Newfile.txt"), and a folder ("~/Test"). 

You can delete files and folders from your computer using the `rm` command.

<div class="alert alert-block alert-danger">
<b>Warning: objects deleted with rm cannot be recovered!</b> They are gone forever, and slight mistakes (like adding a poorly placed space to a file name) could result in the destruction of lots of files. You might even ruin your computer's operating system if you're in the wrong folder. 
</div>

Because of this risk, I recommend that you always inclue the `-i` flag when running `rm`; this calls the "interactive mode", meaning terminal will double-check with you before removing objects.

```
rm -i ~/Desktop/Newfile.txt
```

If you try deleting the folder in the same way it will not work:

```
rm -i ~/Test
```

You will get an error message saying the path `is a directory`. By default folders with files in them cannot be deleted unless you add the `-r` (__recursive__) flag. Recursive mode specifies that you want to delete the folder and all of the files/subfolders within the folder:

````
rm -ir ~/Test
````

## 1.9. Conclusion

That is enough for this lesson. There are many other commands that can be used to navigate around your computer using Terminal. I've put together a cheat sheet below for reference:


## Appendix: Command Line Cheat Sheet
(adapted from https://gist.github.com/poopsplat/7195274#file-gistfile1-textile)CORE 

|command|description|
|-----|-----|
| cd | Home directory |
| cd [folder] | Change directory |
| cd ~ | Home directory, e.g. 'cd ~/folder/' |
| cd / | Root of drive |
| ls | Short listing |
| ls -l | Long listing |
| ls -a | Listing incl. hidden files |
| ls -lh | Long listing with Human readable file sizes |
| ls -R | Entire content of folder recursively |
| sudo [command] | Run command with the security privileges of the superuser (Super User DO) |
| open [file] | Opens a file |
| open . | Opens the directory |
| top | Displays active processes. Press q to quit |
| nano [file] | Opens the Terminal it's editor |
| pico	[file] | Opens the Terminal it's editor |
| q | Exit |
| clear | Clear screen |

### FILE MANAGEMENT

|command|description|
|-----|-----|
| touch [file] | Create new file |
| pwd | Full path to working directory |
| .. | Parent/enclosing directory, e.g. |
| ls -l .. | Long listing of parent directory |
| cd ../../ | Move 2 levels up |
| . | Current folder |
| cat | Concatenate to screen |
| rm [file] | Remove a file, e.g. rm [file] [file] |
| rm -i [file] | Remove with confirmation |
| rm -r [dir] | Remove a directory and contents |
| rm -f [file] | Force removal without confirmation |
| rm -i [file] | Will display prompt before |
| cp [file] [newfile] | Copy file to file |
| cp [file] [dir] | Copy file to directory |
| mv [file] [new filename] | Move/Rename, e.g. mv -v [file] [dir] |

### DIRECTORY MANAGEMENT

|command|description|
|-----|-----|
| mkdir [dir] | Create new directory |
| mkdir -p [dir]/[dir] | Create nested directories |
| rmdir [dir] | Remove directory ( only operates on empty directories ) |
| rm -R [dir] | Remove directory and contents |

### PIPES (Combine multiple commands that generate output)

|command|description|
|-----|-----|
| more | Output content delivered in screensize chunks |
| > [file] | Push output to file, keep in mind it will get overwritten |
| >> [file] | Append output to existing file |
| < | Tell command to read content from a file |

### HELP

|command|description|
|-----|-----|
| [command] -h | Offers help |
| [command] --help | Offers help |
| [command] help | Offers help |
| reset | Resets the terminal display |
| man [command] | Show the help for 'command' |
| whatis [command] | Gives a one-line description of 'command' |
