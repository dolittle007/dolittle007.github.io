
//var b = document.getElementsByTagName("body")[0]
//b.addEventListener("load", passWord, false);

function passWord() {
    var password;
    var pass1="123456";

    password=prompt('Please enter your password to view this page!', ' ');

    if (password==pass1){
       alert('Password Correct! Click OK to enter!');
       window.open("./network.html");
    }
    else{
        window.location="./index.html";
    }

}
