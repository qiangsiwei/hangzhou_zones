# -*- coding: utf-8 -*- 

import math
import glob
import fileinput
import numpy as np
from copy import *
from pylab import *
from scipy import interpolate
import matplotlib.pyplot as plt

# weekend:19,25,26
# weekday:20,21,22,23,24

ranget, rangex, rangey = 24, 200, 200

def euclidean(p1,p2):
	from math import sqrt
	dist = sqrt(sum([pow(i-j,2) for i,j in zip(p1,p2)]))
	return dist

def func1():
	matrix = [[[[] for t in xrange(ranget)] for j in xrange(rangey)] for i in xrange(rangex)]
	for filename in sorted(glob.glob(r"data/*_user_count.txt")):
		print filename
		for line in fileinput.input(filename):
			part = line.strip().split(" ")
			px, py, s, c = int(part[0].split(",")[0]), int(part[0].split(",")[1]), int(part[1]), int(part[2])
			matrix[px][py][s].append(c)
		fileinput.close()
	matrix = [[np.array([np.array(matrix[i][j][t]).sum() for t in xrange(ranget)]).sum()/7 for j in xrange(rangey)] for i in xrange(rangex)]
	out = open("data/all_user_count.txt","w")
	for i in xrange(rangex):
		for j in xrange(rangey):
			out.write(str(i)+" "+str(j)+" "+str(int(matrix[i][j]))+"\n")
	out.close()
	(X, Y), C = meshgrid(np.arange(rangex), np.arange(rangey)), np.array(matrix)
	figure()
	subplot(1,1,1)
	cset1 = pcolormesh(X, Y, C.T, cmap=cm.get_cmap("OrRd"))
	colorbar(cset1)
	axis('off')
	show()
	pass

def func2():
	somlist, sommap = [], {}
	for line in fileinput.input("som/output_si.txt"):
		part = line.strip().replace("[","").replace("]","").split(",")
		px, py = float(part[0]), float(part[1])
		sommap[tuple([px,py])] = 1
	fileinput.close()
	somlist = [[str(i).zfill(4), sommap.keys()[i][0], sommap.keys()[i][1]] for i in xrange(len(sommap.keys()))]
	file = open("som/region.txt","w")
	c = 0
	for line in fileinput.input("data/hz_base.txt"):
		c += 1
		print c
		part = line.strip().split(" ")
		num, lng, lat = part[1]+" "+part[2], float(part[3]), float(part[4])
		if 120.04<=lng<120.24 and 30.16<=lat<30.34:
			minsom, mindist, gx, gy = "", 0, (lng-120.04)/(120.24-120.04)*200, (lat-30.16)/(30.34-30.16)*200
			for som in somlist:
				dist = euclidean(som[1:],[gx, gy])
				if minsom == "" or dist < mindist:
					minsom, mindist = som, dist
			file.write(num+" "+str(gx)+" "+str(gy)+" "+minsom[0]+" "+str(minsom[1])+" "+str(minsom[2])+"\n")
	fileinput.close()
	file.close()
	pass

def func3():
	# # draw grid
	# GR = 50
	# file = open("grid.txt","w")
	# for i in xrange(GR+1):
	# 	lat1, lat2 = 30.16, 30.34
	# 	lng1 = 120.04+i*((120.24-120.04)/GR)
	# 	lng2 = 120.04+i*((120.24-120.04)/GR)
	# 	file.write("map.addOverlay(new BMap.Polyline([new BMap.Point("+str(lng1)+","+str(lat1)+"),new BMap.Point("+str(lng2)+","+str(lat2)+")], {strokeColor:\"red\", strokeWeight:1, strokeOpacity:0.25}));\n")
	# for i in xrange(GR+1):
	# 	lng1, lng2 = 120.04, 120.24
	# 	lat1 = 30.16+i*((30.34-30.16)/GR)
	# 	lat2 = 30.16+i*((30.34-30.16)/GR)
	# 	file.write("map.addOverlay(new BMap.Polyline([new BMap.Point("+str(lng1)+", "+str(lat1)+"),new BMap.Point("+str(lng2)+", "+str(lat2)+")], {strokeColor:\"red\", strokeWeight:1, strokeOpacity:0.25}));\n")
	# file.close()

	# merge spark result
	import fileinput
	seen, num, sommap, assign = set(), 0, {}, {}
	for line in fileinput.input("som/output_si.txt"):
		part = line.strip().replace("[","").replace("]","").split(",")
		if str(part[:2]) not in seen:
			seen.add(str(part[:2]))
			num, px, py = num+1, float(part[0]), float(part[1])
			sommap[num] = {"center":[px,py],"wd_e":[0 for i in xrange(24)],"wd_d":[0 for i in xrange(24)],"we_e":[0 for i in xrange(24)],"we_d":[0 for i in xrange(24)]}
	fileinput.close()
	c = 0
	for line in fileinput.input("data/movement/move_statistic.txt"):
		c += 1
		print c
		part = line.strip().split("\t")[1:]
		for p in part:
			try:
				a, w, i, s = p.split(" ")[0], int(p.split(" ")[1]), int(p.split(" ")[2]), int(p.split(" ")[3])
				if not assign.has_key(tuple([int(a.split(",")[0]),int(a.split(",")[1])])):
					sommin, distmin = 0, 0
					for k,v in sommap.iteritems():
						dist = euclidean(v["center"],[int(a.split(",")[0]),int(a.split(",")[1])])
						if sommin == 0 or dist < distmin:
							sommin, distmin = k, dist
					assign[tuple([int(a.split(",")[0]),int(a.split(",")[1])])] = sommin
				if w == 1 and s == 1:
					sommap[assign[tuple([int(a.split(",")[0]),int(a.split(",")[1])])]]["wd_e"][i] += 1
				if w == 1 and s == 0:
					sommap[assign[tuple([int(a.split(",")[0]),int(a.split(",")[1])])]]["wd_d"][i] += 1
				if w == 0 and s == 1:
					sommap[assign[tuple([int(a.split(",")[0]),int(a.split(",")[1])])]]["we_e"][i] += 1
				if w == 0 and s == 0:
					sommap[assign[tuple([int(a.split(",")[0]),int(a.split(",")[1])])]]["we_d"][i] += 1
			except:
				continue
	fileinput.close()
	file = open("data/movement/move_statistic_merge.txt","w")
	for k,v in sommap.iteritems():
		file.write(str(k)+" "+" ".join([str(i) for i in v["wd_e"]])+" "+" ".join([str(i) for i in v["wd_d"]])+" "+" ".join([str(i) for i in v["we_e"]])+" "+" ".join([str(i) for i in v["we_d"]])+"\n")
	file.close()

	# generate document
	file = open("GibbsLDA++/out.txt","w")
	file.write("368\n")
	c = 0
	for line in fileinput.input("data/movement/move_statistic_merge.txt"):
		part = line.split(" ")
		grid, feature = part[0], [int(i) for i in part[1:]]
		feature = [int(feature[i]/5.0) if 0<=i<24*2 else int(feature[i]/2.0) for i in xrange(len(feature))]
		for i in xrange(24):
			for j in xrange(feature[i]+1):
				file.write(str(i)+",1,E ")
		for i in xrange(24):
			for j in xrange(feature[24+i]+1):
				file.write(str(i)+",1,D ")
		for i in xrange(24):
			for j in xrange(feature[48+i]+1):
				file.write(str(i)+",0,E ")
		for i in xrange(24):
			for j in xrange(feature[72+i]+1):
				file.write(str(i)+",0,D ")
		file.write("\n")
	fileinput.close()
	file.close()

	# ./src/lda -est -alpha 0.5 -beta 0.5 -ntopics 8 -niters 1000 -savestep 2000 -twords 10 -dfile out.txt

	# movement func 1
	seen, somlist = set(), []
	for line in fileinput.input("som/output_si.txt"):
		part = line.strip().replace("[","").replace("]","").split(",")
		if str(part[:2]) not in seen:
			seen.add(str(part[:2]))
			somlist.append([float(part[0])*5, 1000-float(part[1])*5, 0])
	fileinput.close()
	cnt = 0
	for line in fileinput.input("GibbsLDA++/model_movement/model-final.theta"):
		prob = [float(i) for i in line.strip().split(" ")]
		somlist[cnt][2] = prob.index(max(prob))
		cnt += 1
	import pytesselate
	import itertools
	polygons = pytesselate.voronoi(somlist)
	print len(polygons)
	import pydraw
	crs = pydraw.CoordinateSystem([0,0,1000,1000])
	img = pydraw.Image(1000,1000)
	for i,(center,poly) in enumerate(polygons):
		if center:
			if somlist[i][2] == 0:
				img.drawpolygon(poly,fillcolor=(0,0,0),outlinecolor=(0,0,0))
			if somlist[i][2] == 1:
				img.drawpolygon(poly,fillcolor=(255,0,0),outlinecolor=(0,0,0))
			if somlist[i][2] == 2:
				img.drawpolygon(poly,fillcolor=(0,255,0),outlinecolor=(0,0,0))
			if somlist[i][2] == 3:
				img.drawpolygon(poly,fillcolor=(0,0,255),outlinecolor=(0,0,0))
			if somlist[i][2] == 4:
				img.drawpolygon(poly,fillcolor=(0,255,255),outlinecolor=(0,0,0))
			if somlist[i][2] == 5:
				img.drawpolygon(poly,fillcolor=(255,0,255),outlinecolor=(0,0,0))
			if somlist[i][2] == 6:
				img.drawpolygon(poly,fillcolor=(255,255,0),outlinecolor=(0,0,0))
			if somlist[i][2] == 7:
				img.drawpolygon(poly,fillcolor=(255,255,255),outlinecolor=(0,0,0))
			# img.drawsquare(*center[:2], fillsize=2, fillcolor=(255,0,0))
	img.save("test.png")
	img.view()
	# y0list, y1list, y2list, y3list = [], [], [], []
	# # for line in fileinput.input("GibbsLDA++/model-final.phi"):
	# for line in fileinput.input("LDA/model-final.phi"):
	# 	prob = [float(i) for i in line.strip().split(" ")]
	# 	y0, y1, y2, y3 = prob[0:24], prob[24:48], prob[48:72], prob[72:96]
	# 	y0list.append(y0)
	# 	y1list.append(y1)
	# 	y2list.append(y2)
	# 	y3list.append(y3)
	# fig = plt.figure()
	# ax1 = fig.add_subplot(241)
	# ax1.plot([i for i in xrange(24)], y0list[0], c="red")
	# ax1.plot([i for i in xrange(24)], y1list[0], c="yellow")
	# ax1.plot([i for i in xrange(24)], y2list[0], c="green")
	# ax1.plot([i for i in xrange(24)], y3list[0], c="blue")
	# ax1.set_xlim(0,23)
	# ax1.set_ylim(0,0.15)
	# ax2 = fig.add_subplot(242)
	# ax2.plot([i for i in xrange(24)], y0list[1], c="red")
	# ax2.plot([i for i in xrange(24)], y1list[1], c="yellow")
	# ax2.plot([i for i in xrange(24)], y2list[1], c="green")
	# ax2.plot([i for i in xrange(24)], y3list[1], c="blue")
	# ax2.set_xlim(0,23)
	# ax2.set_ylim(0,0.15)
	# ax3 = fig.add_subplot(243)
	# ax3.plot([i for i in xrange(24)], y0list[2], c="red")
	# ax3.plot([i for i in xrange(24)], y1list[2], c="yellow")
	# ax3.plot([i for i in xrange(24)], y2list[2], c="green")
	# ax3.plot([i for i in xrange(24)], y3list[2], c="blue")
	# ax3.set_xlim(0,23)
	# ax3.set_ylim(0,0.15)
	# ax4 = fig.add_subplot(244)
	# ax4.plot([i for i in xrange(24)], y0list[3], c="red")
	# ax4.plot([i for i in xrange(24)], y1list[3], c="yellow")
	# ax4.plot([i for i in xrange(24)], y2list[3], c="green")
	# ax4.plot([i for i in xrange(24)], y3list[3], c="blue")
	# ax4.set_xlim(0,23)
	# ax4.set_ylim(0,0.15)
	# ax5 = fig.add_subplot(245)
	# ax5.plot([i for i in xrange(24)], y0list[4], c="red")
	# ax5.plot([i for i in xrange(24)], y1list[4], c="yellow")
	# ax5.plot([i for i in xrange(24)], y2list[4], c="green")
	# ax5.plot([i for i in xrange(24)], y3list[4], c="blue")
	# ax5.set_xlim(0,23)
	# ax5.set_ylim(0,0.15)
	# ax6 = fig.add_subplot(246)
	# ax6.plot([i for i in xrange(24)], y0list[5], c="red")
	# ax6.plot([i for i in xrange(24)], y1list[5], c="yellow")
	# ax6.plot([i for i in xrange(24)], y2list[5], c="green")
	# ax6.plot([i for i in xrange(24)], y3list[5], c="blue")
	# ax6.set_xlim(0,23)
	# ax6.set_ylim(0,0.15)
	# ax7 = fig.add_subplot(247)
	# ax7.plot([i for i in xrange(24)], y0list[6], c="red")
	# ax7.plot([i for i in xrange(24)], y1list[6], c="yellow")
	# ax7.plot([i for i in xrange(24)], y2list[6], c="green")
	# ax7.plot([i for i in xrange(24)], y3list[6], c="blue")
	# ax7.set_xlim(0,23)
	# ax7.set_ylim(0,0.15)
	# ax8 = fig.add_subplot(248)
	# ax8.plot([i for i in xrange(24)], y0list[7], c="red")
	# ax8.plot([i for i in xrange(24)], y1list[7], c="yellow")
	# ax8.plot([i for i in xrange(24)], y2list[7], c="green")
	# ax8.plot([i for i in xrange(24)], y3list[7], c="blue")
	# ax8.set_xlim(0,23)
	# ax8.set_ylim(0,0.15)
	# plt.show()

	# movement func 2
	seen, somlist = set(), []
	for line in fileinput.input("som/output_si.txt"):
		part = line.strip().replace("[","").replace("]","").split(",")
		if str(part[:2]) not in seen:
			seen.add(str(part[:2]))
			somlist.append([float(part[0])*5, 1000-float(part[1])*5, [0]*8])
	fileinput.close()
	print len(somlist)
	cnt = 0
	for line in fileinput.input("GibbsLDA++/model-final.theta"):
		somlist[cnt][2] = [float(i) for i in line.strip().split(" ")]
		cnt += 1
	import pytesselate
	import itertools
	polygons = pytesselate.voronoi(somlist)
	import pydraw
	crs = pydraw.CoordinateSystem([0,0,1000,1000])
	img = pydraw.Image(1000,1000)
	c = 7
	for i,(center,poly) in enumerate(polygons):
		if center:
			img.drawpolygon(poly,fillcolor=(255*somlist[i][2][c],255*somlist[i][2][c],255*somlist[i][2][c]),outlinecolor=(0,0,0))
	img.save("figure/class_"+str(c)+".png")
	# img.view()

	# # poi func1
	# tagmap = {}
	# for line in fileinput.input("data/poi/poi.txt"):
	# 	part = line.strip().split("\t")
	# 	tag, lng, lat = part[0], float(part[1]), float(part[2])
	# 	if 120.04<=lng<=120.24 and 30.16<=lat<=30.34:
	# 		tagmap[tag] = 1 if not tagmap.has_key(tag) else tagmap[tag]+1
	# fileinput.close()
	# taglist = []
	# for k,v in tagmap.iteritems():
	# 	taglist.append({"k":k,"v":v})
	# taglist = sorted(taglist, key=lambda x:x["v"])
	# tag, count = [], []
	# for item in taglist:
	# 	if item["k"] in ["公司","购物","美食","房产","休闲","政府","丽人","金融","医疗","宾馆","汽车","教育","培训","旅游","运动"]:
	# 		item["k"] = "company" if item["k"] == "公司"\
	# 		else "shopping" if item["k"] == "购物"\
	# 		else "food" if item["k"] == "美食"\
	# 		else "housing" if item["k"] == "房产"\
	# 		else "leisure" if item["k"] == "休闲"\
	# 		else "government" if item["k"] == "政府"\
	# 		else "beauty" if item["k"] == "丽人"\
	# 		else "finance" if item["k"] == "金融"\
	# 		else "medical" if item["k"] == "医疗"\
	# 		else "hotel" if item["k"] == "宾馆"\
	# 		else "car" if item["k"] == "汽车"\
	# 		else "education" if item["k"] == "教育"\
	# 		else "training" if item["k"] == "培训"\
	# 		else "traveling" if item["k"] == "旅游"\
	# 		else "physical" if item["k"] == "运动"\
	# 		else ""
	# 		tag.append(unicode(item["k"],"utf-8"))
	# 		count.append(item["v"])
	# for i in xrange(len(count)):
	# 	print round(float(count[i])/sum(count),4)
	# import matplotlib.pyplot as plt
	# import numpy as np
	# y_pos = np.arange(len(tag))
	# plt.barh(y_pos, count, align='center', alpha=0.4)
	# plt.yticks(y_pos, tag)
	# plt.xlabel('count')
	# plt.title('POI statistics')
	# plt.show()

	# # poi func2
	# seen, num, sommap, assign = set(), 0, {}, {}
	# for line in fileinput.input("som/output_si.txt"):
	# 	part = line.strip().replace("[","").replace("]","").split(",")
	# 	if str(part[:2]) not in seen:
	# 		seen.add(str(part[:2]))
	# 		px, py = float(part[0]), float(part[1])
	# 		sommap[num] = {"center":[px,py],"poi":[0]*15}
	# 		num += 1
	# fileinput.close()
	# for line in fileinput.input("data/poi/poi.txt"):
	# 	part = line.strip().split("\t")
	# 	tag, lng, lat = part[0], float(part[1]), float(part[2])
	# 	if tag in ["公司","购物","美食","房产","休闲","政府","丽人","金融","医疗","宾馆","汽车","教育","培训","旅游","运动"] and 120.04<=lng<120.24 and 30.16<=lat<30.34:
	# 		lng, lat = int((lng-120.04)/(120.24-120.04)*200), int((lat-30.16)/(30.34-30.16)*200)
	# 		tag = 0 if tag == "公司"\
	# 		else  1 if tag == "购物"\
	# 		else  2 if tag == "美食"\
	# 		else  3 if tag == "房产"\
	# 		else  4 if tag == "休闲"\
	# 		else  5 if tag == "政府"\
	# 		else  6 if tag == "丽人"\
	# 		else  7 if tag == "金融"\
	# 		else  8 if tag == "医疗"\
	# 		else  9 if tag == "宾馆"\
	# 		else 10 if tag == "汽车"\
	# 		else 11 if tag == "教育"\
	# 		else 12 if tag == "培训"\
	# 		else 13 if tag == "旅游"\
	# 		else 14 if tag == "运动"\
	# 		else -1
	# 		if not assign.has_key(tuple([lng, lat])):
	# 			sommin, distmin = 0, 0
	# 			for k,v in sommap.iteritems():
	# 				dist = euclidean(v["center"],[lng, lat])
	# 				if sommin == 0 or dist < distmin:
	# 					sommin, distmin = k, dist
	# 			assign[tuple([lng, lat])] = sommin
	# 		sommap[assign[tuple([lng, lat])]]["poi"][tag] += 1
	# fileinput.close()
	# classlist = [[0]*15 for i in xrange(8)]
	# num = 0
	# for line in fileinput.input("GibbsLDA++/model-final.theta"):
	# 	prob = [float(i) for i in line.strip().split(" ")]
	# 	for c in xrange(8):
	# 		classlist[c] = [i+j for i,j in zip(classlist[c], [k*prob[c] for k in sommap[num]["poi"]])]
	# 	num += 1
	# fileinput.close()
	# print classlist
	# classlist = [[float(classlist[i][j])/sum(classlist[i]) for j in xrange(15)] for i in xrange(8)]
	# avg = [0.2304,0.2043,0.1472,0.1332,0.0504,0.0393,0.0391,0.0321,0.0263,0.0242,0.0223,0.0185,0.0128,0.0100,0.0099]
	# classlist = [[classlist[i][j]-avg[j] for j in xrange(15)] for i in xrange(8)]
	# import matplotlib.pyplot as plt
	# import numpy as np
	# fig = plt.figure()
	# ax1 = fig.add_subplot(241)
	# tag = [str(i) for i in xrange(15)]
	# y_pos = np.arange(len(tag))
	# plt.bar(y_pos, classlist[0], align='center', alpha=0.4)
	# plt.xticks(y_pos, tag)
	# ax1.set_ylim(-0.1,0.1)
	# ax2 = fig.add_subplot(242)
	# tag = [str(i) for i in xrange(15)]
	# y_pos = np.arange(len(tag))
	# plt.bar(y_pos, classlist[1], align='center', alpha=0.4)
	# plt.xticks(y_pos, tag)
	# ax2.set_ylim(-0.1,0.1)
	# ax3 = fig.add_subplot(243)
	# tag = [str(i) for i in xrange(15)]
	# y_pos = np.arange(len(tag))
	# plt.bar(y_pos, classlist[2], align='center', alpha=0.4)
	# plt.xticks(y_pos, tag)
	# ax3.set_ylim(-0.1,0.1)
	# ax4 = fig.add_subplot(244)
	# tag = [str(i) for i in xrange(15)]
	# y_pos = np.arange(len(tag))
	# plt.bar(y_pos, classlist[3], align='center', alpha=0.4)
	# plt.xticks(y_pos, tag)
	# ax4.set_ylim(-0.1,0.1)
	# ax5 = fig.add_subplot(245)
	# tag = [str(i) for i in xrange(15)]
	# y_pos = np.arange(len(tag))
	# plt.bar(y_pos, classlist[4], align='center', alpha=0.4)
	# plt.xticks(y_pos, tag)
	# ax5.set_ylim(-0.1,0.1)
	# ax6 = fig.add_subplot(246)
	# tag = [str(i) for i in xrange(15)]
	# y_pos = np.arange(len(tag))
	# plt.bar(y_pos, classlist[5], align='center', alpha=0.4)
	# plt.xticks(y_pos, tag)
	# ax6.set_ylim(-0.1,0.1)
	# ax7 = fig.add_subplot(247)
	# tag = [str(i) for i in xrange(15)]
	# y_pos = np.arange(len(tag))
	# plt.bar(y_pos, classlist[6], align='center', alpha=0.4)
	# plt.xticks(y_pos, tag)
	# ax7.set_ylim(-0.1,0.1)
	# ax8 = fig.add_subplot(248)
	# tag = [str(i) for i in xrange(15)]
	# y_pos = np.arange(len(tag))
	# plt.bar(y_pos, classlist[7], align='center', alpha=0.4)
	# plt.xticks(y_pos, tag)
	# ax8.set_ylim(-0.1,0.1)
	# plt.show()

	# generate document poi
	seen, num, sommap, assign = set(), 0, {}, {}
	for line in fileinput.input("som/output_si.txt"):
		part = line.strip().replace("[","").replace("]","").split(",")
		if str(part[:2]) not in seen:
			seen.add(str(part[:2]))
			px, py = float(part[0]), float(part[1])
			sommap[num] = {"center":[px,py],"poi":[0]*15}
			num += 1
	fileinput.close()
	count = 0
	for line in fileinput.input("data/poi/poi.txt"):
		count += 1
		print count
		part = line.strip().split("\t")
		tag, lng, lat = part[0], float(part[1]), float(part[2])
		if tag in ["公司","购物","美食","房产","休闲","政府","丽人","金融","医疗","宾馆","汽车","教育","培训","旅游","运动"] and 120.04<=lng<120.24 and 30.16<=lat<30.34:
			lng, lat = int((lng-120.04)/(120.24-120.04)*200), int((lat-30.16)/(30.34-30.16)*200)
			tag = 0 if tag == "公司"\
			else  1 if tag == "购物"\
			else  2 if tag == "美食"\
			else  3 if tag == "房产"\
			else  4 if tag == "休闲"\
			else  5 if tag == "政府"\
			else  6 if tag == "丽人"\
			else  7 if tag == "金融"\
			else  8 if tag == "医疗"\
			else  9 if tag == "宾馆"\
			else 10 if tag == "汽车"\
			else 11 if tag == "教育"\
			else 12 if tag == "培训"\
			else 13 if tag == "旅游"\
			else 14 if tag == "运动"\
			else -1
			if not assign.has_key(tuple([lng, lat])):
				sommin, distmin = 0, 0
				for k,v in sommap.iteritems():
					dist = euclidean(v["center"],[lng, lat])
					if sommin == 0 or dist < distmin:
						sommin, distmin = k, dist
				assign[tuple([lng, lat])] = sommin
			sommap[assign[tuple([lng, lat])]]["poi"][tag] += 1
	fileinput.close()
	file = open("GibbsLDA++/out_poi.txt","w")
	file.write("368\n")
	for k,v in sommap.iteritems():
		file.write("-1 "+" ".join([" ".join([str(x) for y in xrange(v["poi"][x])]) for x in xrange(15)])+"\n")
	file.close()

	# ./src/lda -est -alpha 0.5 -beta 0.5 -ntopics 8 -niters 1000 -savestep 2000 -twords 10 -dfile out_poi.txt

	# movement func 1
	seen, somlist = set(), []
	for line in fileinput.input("som/output_si.txt"):
		part = line.strip().replace("[","").replace("]","").split(",")
		if str(part[:2]) not in seen:
			seen.add(str(part[:2]))
			somlist.append([float(part[0])*5, 1000-float(part[1])*5, 0])
	fileinput.close()
	cnt = 0
	for line in fileinput.input("GibbsLDA++/model-final.theta"):
		prob = [float(i) for i in line.strip().split(" ")]
		somlist[cnt][2] = prob.index(max(prob))
		cnt += 1
	import pytesselate
	import itertools
	polygons = pytesselate.voronoi(somlist)
	print len(polygons)
	import pydraw
	crs = pydraw.CoordinateSystem([0,0,1000,1000])
	img = pydraw.Image(1000,1000)
	for i,(center,poly) in enumerate(polygons):
		if center:
			if somlist[i][2] == 0:
				img.drawpolygon(poly,fillcolor=(0,0,0),outlinecolor=(0,0,0))
			if somlist[i][2] == 1:
				img.drawpolygon(poly,fillcolor=(255,0,0),outlinecolor=(0,0,0))
			if somlist[i][2] == 2:
				img.drawpolygon(poly,fillcolor=(0,255,0),outlinecolor=(0,0,0))
			if somlist[i][2] == 3:
				img.drawpolygon(poly,fillcolor=(0,0,255),outlinecolor=(0,0,0))
			if somlist[i][2] == 4:
				img.drawpolygon(poly,fillcolor=(0,255,255),outlinecolor=(0,0,0))
			if somlist[i][2] == 5:
				img.drawpolygon(poly,fillcolor=(255,0,255),outlinecolor=(0,0,0))
			if somlist[i][2] == 6:
				img.drawpolygon(poly,fillcolor=(255,255,0),outlinecolor=(0,0,0))
			if somlist[i][2] == 7:
				img.drawpolygon(poly,fillcolor=(255,255,255),outlinecolor=(0,0,0))
			# img.drawsquare(*center[:2], fillsize=2, fillcolor=(255,0,0))
	img.save("test.png")
	img.view()
	# ylist = []
	# for line in fileinput.input("GibbsLDA++/model-final.phi"):
	# 	prob = [float(i) for i in line.strip().split(" ")]
	# 	ylist.append(prob)
	# print len(ylist[0])
	# fig = plt.figure()
	# ax1 = fig.add_subplot(241)
	# ax1.plot([i for i in xrange(16)], ylist[0], c="red")
	# ax1.set_xlim(0,15)
	# ax1.set_ylim(0,1.0)
	# ax2 = fig.add_subplot(242)
	# ax2.plot([i for i in xrange(16)], ylist[1], c="red")
	# ax2.set_xlim(0,15)
	# ax2.set_ylim(0,1.0)
	# ax3 = fig.add_subplot(243)
	# ax3.plot([i for i in xrange(16)], ylist[2], c="red")
	# ax3.set_xlim(0,15)
	# ax3.set_ylim(0,1.0)
	# ax4 = fig.add_subplot(244)
	# ax4.plot([i for i in xrange(16)], ylist[3], c="red")
	# ax4.set_xlim(0,15)
	# ax4.set_ylim(0,1.0)
	# ax5 = fig.add_subplot(245)
	# ax5.plot([i for i in xrange(16)], ylist[4], c="red")
	# ax5.set_xlim(0,15)
	# ax5.set_ylim(0,1.0)
	# ax6 = fig.add_subplot(246)
	# ax6.plot([i for i in xrange(16)], ylist[5], c="red")
	# ax6.set_xlim(0,15)
	# ax6.set_ylim(0,1.0)
	# ax7 = fig.add_subplot(247)
	# ax7.plot([i for i in xrange(16)], ylist[6], c="red")
	# ax7.set_xlim(0,15)
	# ax7.set_ylim(0,1.0)
	# ax8 = fig.add_subplot(248)
	# ax8.plot([i for i in xrange(16)], ylist[7], c="red")
	# ax8.set_xlim(0,15)
	# ax8.set_ylim(0,1.0)
	# plt.show()

	# movement func 2
	seen, somlist = set(), []
	for line in fileinput.input("som/output_si.txt"):
		part = line.strip().replace("[","").replace("]","").split(",")
		if str(part[:2]) not in seen:
			seen.add(str(part[:2]))
			somlist.append([float(part[0])*5, 1000-float(part[1])*5, [0]*8])
	fileinput.close()
	print len(somlist)
	cnt = 0
	for line in fileinput.input("GibbsLDA++/model-final.theta"):
		somlist[cnt][2] = [float(i) for i in line.strip().split(" ")]
		cnt += 1
	import pytesselate
	import itertools
	polygons = pytesselate.voronoi(somlist)
	import pydraw
	crs = pydraw.CoordinateSystem([0,0,1000,1000])
	img = pydraw.Image(1000,1000)
	c = 6
	for i,(center,poly) in enumerate(polygons):
		if center:
			img.drawpolygon(poly,fillcolor=(255*somlist[i][2][c],255*somlist[i][2][c],255*somlist[i][2][c]),outlinecolor=(0,0,0))
	img.save("figure/class_"+str(c)+".png")
	# img.view()

if __name__ == "__main__":
	func1()
	func2()
	func3()
