load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/cd_string.ncl"

if (.not. isdefined("file_prefix")) then
	file_prefix = "rh_test"
end if

if (isdefined("time_range")) then

else
  file_paths = systemfunc("ls " + file_prefix + "*.nc")
end if

if (.not. isdefined("wks_type")) then
	wks_type = "pdf"
end if
wks_type@wkWidth = 500
wks_type@wkHeight = 500

wks = gsn_open_wks(wks_type, file_prefix)

gsn_define_colormap(wks, "BlAqGrYeOrRe")

res = True
res@gsnDraw = False
res@gsnFrame = False

res_gd = res
res_gd@cnLineLabelsOn = False
res_gd@cnLineThicknessF = 1.0
res_gd@cnLevelSelectionMode = "ManualLevels"
res_gd@cnMinLevelValF = 7.5e4
res_gd@cnMaxLevelValF = 10.5e4
res_gd@cnLevelSpacingF = 0.2e4
res_gd@gsnPolar = "NH"

res_text = True
res_text@txFontHeightF = 0.015

res_anno = True
res_anno@amParallelPosF = -0.5
res_anno@amOrthogonalPosF = 0.6

do k = 0, dimsizes(file_paths) - 1
  print("[Notice]: Plot " + file_paths(k))
  f = addfile(file_paths(k), "r")

	params = getvaratts(f)
	text = ""
	do i = 0, dimsizes(params) - 1
		if (params(i) .eq. "dataset" .or. params(i) .eq. "author") then
			continue
		end if
		text = text + "~C~" + params(i) + ": " + f@$params(i)$
	end do
  
	time = f->time
  time@units = "hours since 0000-01-01 00:00:00"
  res_gd@gsnRightString = cd_string(f->time, "%J %H:%M")
  plot = gsn_csm_contour_map(wks, f->gd(0,:,:), res_gd)
	anno = gsn_add_annotation(plot, gsn_create_text(wks, text, res_text), res_anno)
  draw(plot)
  frame(wks)
end do
