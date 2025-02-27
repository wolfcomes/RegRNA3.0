function Window() {
  var Win = null;

  var nArgs = Window.arguments.length;
  var Url = Window.arguments[0];
  var Target = Window.arguments[1];
  var Width = Window.arguments[2];
  var Height = Window.arguments[3];

	if (nArgs < 4) {
		Width = SGL_JS_WINWIDTH;
		Height = SGL_JS_WINHEIGHT;
	}

	Win = window.open('', Target,
		'toobar=no,location=no,directories=no,status=yes,scrollbars=yes,resizable=yes,menubar=no,copyhistory=no,' +
		'width=' + Width + ',height=' + Height);

	Win.location.href = Url;

	Win.focus();
}

// ------------------------------------------------------------------------------------------------
