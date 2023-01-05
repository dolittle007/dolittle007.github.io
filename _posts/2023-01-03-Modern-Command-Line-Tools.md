---
layout: post
title: "Modern commnad line tools powered by Rust"
date: 2023-01-03
category: how-to
tags: [command line tools, Rust]
---

An introduction to modern command line tools powered by Rust.

<!--more-->

Rust is taking over the terminal. Rust is a general-purpose programming language that is blazing fast and memory safe. It is the fastest-growing and most loved programming language in the world. It is used to build everything from operating systems to web servers to command-line tools. Recently there has been a surge of command line tools and utilities written in Rust, and many of them are intended to replace standard Unix commands. They are faster, more user-friendly, and have more features than their standard Unix counterparts. In this post, I will cover some of the best Rust command line tools I have used for a while. You can also use these to supercharge your terminal.

### Alacritty

Let us start with the terminal itself. [Alacritty](https://github.com/alacritty/alacritty) is a cross-platform modern terminal emulator with sensible defaults. It is GPU accelerated, super fast, and highly configurable. You can use it on Linux, macOS, and Windows.


### Installation
```bash
# Arch Linux
yay -S alacritty
# Fedora/CentOS
dnf copr enable atim/alacritty
dnf install alacritty
# Debian/Ubuntu
add-apt-repository ppa:aslatter/ppa
apt install alacritty
# macOS Homebrew
brew install --cask alacritty
# Windows Scoop
scoop bucket add extras
scoop install alacritty
# Cargo on any
cargo install alacritty
```

### Starship

[Starship](https://starship.rs/) is the best terminal prompt I have ever used. Forget Oh My Zsh and stuff like that. Starship is fast, highly customizable, and has a great default theme and settings. I didn’t even change most of the default settings, as things were perfect as it is. Starship works on shells like zsh, fish, and bash and can also work alongside other prompts like Oh My Zsh, in case you still want to use Oh My Zsh for other plugins like autosuggestions and so on. Starship works best with a Nerd Font as it can show icons and ligatures based on context. I used Oh My Zsh for many years with the powerlevel10k theme, but the prompt was a bit slow. Starship is blazing fast with more features and an excellent UX

### Installation
```bash
# Arch Linux
yay -S starship
# Fedora/CentOS
dnf install starship
# Debian/Ubuntu
curl -sS https://starship.rs/install.sh | sh
# macOS/Linux Homebrew
brew install starship
# macOS MacPorts
port install starship
# Windows Scoop
scoop install starship
# Cargo
cargo install starship --locked
```

### bat

[bat](https://github.com/sharkdp/bat) is one of my favorite tools from this list. It’s a replacement for cat, and once you have used bat, you will never go back.

### Installation
```bash
# Arch Linux
yay -S bat
# Fedora/CentOS
dnf install bat
# Debian/Ubuntu
apt install bat
# macOS/Linux Homebrew
brew install bat
# macOS MacPorts
port install bat
# Windows Scoop
scoop install bat
# Cargo
cargo install bat --locked
```

### rip

rip is an improved version of the rm command. It is faster, safer, and user-friendly. rip sends deleted files to a temp location so they can be recovered using rip -u. I really like the simplicity and the revert feature, as I don’t have to worry about accidentally deleting something using rm.

### Installation

```bash
# Arch Linux
yay -S rm-improved
# Fedora/CentOS/Debian/Ubuntu
# Install from binary or build locally using Cargo
# macOS Homebrew
brew install rm-improved
# Cargo
cargo install rm-improved
```

### zoxide

zoxide is a smarter cd replacement.

### Installation

```bash
# Arch Linux
yay -S zoxide
# Fedora/CentOS
dnf install zoxide
# Debian/Ubuntu
apt install zoxide
# macOS/Linux Homebrew
brew install zoxide
# macOS MacPorts
port install zoxide
# Windows Scoop
scoop install zoxide
# Cargo
cargo install zoxide --locked
```

Once installed, you must add the following to your shell config file.

```bash
# bash (~/.bashrc)
eval "$(zoxide init bash)"
# zsh (~/.zshrc)
eval "$(zoxide init zsh)"
# fish (~/.config/fish/config.fish)
zoxide init fish | source

# Alias cd to z
alias cd='z'
```
### LSD and exa

Both [LSD](https://github.com/Peltoche/lsd) and [exa](https://github.com/ogham/exa) are replacements for the ls command. They both look gorgeous with nice colors and icons and have features like headers, sorting, tree views, and so on. Exa is a bit faster than LSD for tree views and can show the Git status of files and folders.

### exa installation
```bash
# Arch Linux
yay -S exa
# Fedora/CentOS
dnf install exa
# Debian/Ubuntu
apt install exa
# macOS Homebrew
brew install exa
# Cargo
cargo install exa

# Alias ls to exa
alias ls='exa --git --icons --color=always --group-directories-first'

```
### LSD installation
```bash
# Arch Linux
yay -S lsd
# Fedora/CentOS
dnf install lsd
# Debian/Ubuntu
dpkg -i lsd_0.23.1_amd64.deb # get .deb file from https://github.com/Peltoche/lsd/releases
# macOS Homebrew
brew install lsd
# macOS MacPorts
port install lsd
# Windows Scoop
scop install lsd
# Cargo
cargo install lsd

# Alias ls to lsd
alias ls='lsd --header --color=always --group-directories-first'
```

### ripgrep
[ripgrep (rg)](https://github.com/BurntSushi/ripgrep) is a line-oriented search tool that recursively searches your current directory for a regex pattern. It is faster than grep and has many features like compressed files search, colorized output, smart case, file type filtering, multi-threading, and so on.

### Installation

```bash
# Arch Linux
yay -S ripgrep
# Fedora/CentOS
dnf install ripgrep
# Debian/Ubuntu
apt-get install ripgrep
# macOS/Linux Homebrew
brew install ripgrep
# macOS MacPorts
port install ripgrep
# Windows Scoop
scoop install ripgrep
# Cargo
cargo install ripgrep
```

### fd
[fd](https://github.com/sharkdp/fd) is a simpler alternative to __find__.  It is extremely fast due to parallel traversing and shows a modern colorized output and supports patterns and regex, parallel commands, smart case.

### Installation
```bash
# Arch Linux
yay -S fd
# Fedora/CentOS
dnf install fd-find
# Debian/Ubuntu
apt install fd-find
# macOS Homebrew
brew install fd
# macOS MacPorts
port install fd
# Windows Scoop
scoop install fd
# Cargo
cargo install fd-find
```

### dust
[Dust](https://github.com/bootandy/dust) is an alternative for the du command. It is fast and has a better UX with nice visualization for disk usage.

### Installation
```bash
# Arch Linux
yay -S dust
# Fedora/CentOS
# Install binary from https://github.com/bootandy/dust/releases
# Debian/Ubuntu
deb-get install du-dust
# macOS Homebrew
brew install dust
# macOS MacPorts
port install dust
# Windows Scoop
scoop install dust
# Cargo
cargo install du-dust
```

### Hyperfine

Hyperfine is a Rust-powered, __time__ alternative, command-line benchmarking tool. 

### Installation
```bash
# Arch Linux
yay -S hyperfine
# Fedora/CentOS
# Install binary from https://github.com/bootandy/dust/releases
# Debian/Ubuntu
deb-get install hyperfine
# macOS Homebrew
brew install hyperfine
# macOS MacPorts
port install hyperfine
# Windows Scoop
scoop install hyperfine
# Cargo
cargo install hyperfine
```

### sd

[sd](https://github.com/chmln/sd) is a find-and-replace CLI, and you can use it as a replacement for __sed__ and __awk__.

### Installation
```bash
# Arch Linux
yay -S sd
# Fedora/CentOS
dnf install sd
# Debian/Ubuntu
# Install binary from the release page
# macOS Homebrew
brew install sd
# Windows Scoop
choco install sd-cli
# Cargo
cargo install sd
```


### bottom
[bottom](https://github.com/ClementTsang/bottom) is a __top__ replacement with a nice terminal UI. It’s quite feature-rich and customizable.

### Installation
```bash
# Arch Linux
yay -S bottom
# Fedora/CentOS
dnf copr enable atim/bottom -y
dnf install bottom
# Debian/Ubuntu
dpkg -i bottom_0.6.8_amd64.deb
# macOS Homebrew
brew install bottom
# macOS MacPorts
port install bottom
# Windows Scoop
scoop install bottom
# Cargo
cargo install bottom --locked

```

### A Bonus Tip

When you use another computer, server, or system, you will be using Linux commands. It is a good idea to keep using Linux commands even if you are using Rust-powered alternatives.

I have the following aliases in my .zshrc and .bashrc.

```bash
alias ls='exa'
alias x='exa'
alias find ='fd'
alias grep='rg'
alias du='dust'
alias cat='bat'
alias time='hyperfine'
alias sed='sd'
alias top='btm'
alias htop='btm'
```

You can find a list of other Rust CLI tools [here](https://gist.github.com/sts10/daadbc2f403bdffad1b6d33aff016c0a).


### References
* [Rust Easy!](https://deepu.tech/rust-terminal-tools-linux-mac-windows-fish-zsh/)
* [A curated list of command-line utilities written in Rust](https://gist.github.com/sts10/daadbc2f403bdffad1b6d33aff016c0a)
* [A collection of notable Rust blog posts](https://gist.github.com/brson/a324c83a6af6a8a78dfaa9d33eb9b48e)

