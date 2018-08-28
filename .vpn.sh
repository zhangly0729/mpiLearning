#!/usr/bin/expect -f
# 设置ssh连接的用户名
set user wsh18
# 设置ssh连接的host地址
set host 41.0.0.188
# 设置ssh连接的port端口号
set port 22
# 设置ssh连接的登录密码
set password Awusihai18 
# 设置ssh连接的超时时间
set timeout -1

spawn ssh $user@$host -p $port
#expect "(yes/no)"
#send "yes\r";
expect "*password:"
# 提交密码
send "$password\r"
send "cd zsh"
send "cd online1"

#控制权移交
interact

