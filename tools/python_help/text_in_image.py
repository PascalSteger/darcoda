import Image, ImageDraw, ImageFont
from os import chdir, path

def txt2img(fi,text,bg="#ffffff",fg="#000000",font="FreeSans.ttf",FontSize=14):
    font_dir = "/usr/share/fonts/truetype/freefont/"
    img_name = fi+".jpg"
    font_size = FontSize
    fnt = ImageFont.truetype(font_dir+font, font_size)
    lineWidth = 20
    img = Image.open(fi)
    # make an entirely black image
    imgbg = Image.new('RGBA', img.size, "#000000") 
    mask = Image.new('L',img.size,"#000000")       # make a mask that masks out all
    draw = ImageDraw.Draw(img)                     # setup to draw on the main image
    drawmask = ImageDraw.Draw(mask)                # setup to draw on the mask
    drawmask.line((0, lineWidth/2, img.size[0],lineWidth/2),
                  fill="#999999", width=10)        # draw a line on the mask to allow some bg through
    img.paste(imgbg, mask=mask)                    # put the (somewhat) transparent bg on the main
    draw.text((10,0), text, font=fnt, fill=bg)      # add some text to the main
    del draw
    img.save(img_name,"JPEG",quality=100)
    
txt2img("pic.png","This is a really weird image")
