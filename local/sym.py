from shutil import copy
from time import sleep
from hashlib import md5

while True:
    try:
        with open(r"C:/Users/ZH/Desktop/汇报/Janus WSSe单层：一种出色的光解水催化剂.docx", "rb") as f1:
            with open(r"C:/Users/ZH/OneDrive - stu.suda.edu.cn/Janus WSSe单层：一种出色的光解水催化剂.docx", "rb") as f2:
                content1 = f1.read()
                content2 = f2.read()
                # print(md5(content1).hexdigest())
                # print(md5(content2).hexdigest())
            if md5(content1).hexdigest() == md5(content2).hexdigest():
                print("未更改进入休眠！")
                f1.close()
                f2.close()
                sleep(2)
            else:
                copy(r"C:/Users/ZH/Desktop/汇报/Janus WSSe单层：一种出色的光解水催化剂.docx", r"C:/Users/ZH/OneDrive - stu.suda.edu.cn")
                print("同步成功！")
                continue
    except PermissionError:
        print("权限出现错误,请勿打开‘C:/Users/ZH/Desktop/汇报/Janus WSSe单层：一种出色的光解水催化剂.docx文件’")
        sleep(2)
        continue
