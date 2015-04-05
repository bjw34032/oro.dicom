context("Does testthat even work?")

test_that("Convert", {
  expect_equal(dec2base(23, 2), "10111")
  expect_equal(dec2base(2, 2), "10")
})


abdo <- system.file("dcm/Abdo.dcm", package="oro.dicom")

test_that("Reading DICOM file Abdo.dcm", {
  expect_is(readDICOMFile(abdo), "list")
  expect_is(readDICOMFile(abdo)$hdr, "data.frame")
  expect_is(readDICOMFile(abdo)$img, "matrix")
})

spine <- system.file("dcm/Spine1.dcm", package="oro.dicom")

test_that("Reading DICOM file Spine1.dcm", {
  expect_is(readDICOMFile(spine), "list")
  expect_is(readDICOMFile(spine)$hdr, "data.frame")
  expect_is(readDICOMFile(spine)$img, "matrix")
})

sphere <- system.file("sphere3", package="oro.dicom")

test_that("Reading DICOM files sphere3", {
  expect_is(readDICOM(sphere), "list")
  expect_is(readDICOM(sphere)[[1]]$hdr, "vector")
  expect_is(readDICOM(sphere)[[1]]$img, "matrix")
})
