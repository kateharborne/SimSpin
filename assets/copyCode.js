// This assumes that you're using Rouge; if not, update the selector
// const codeBlocks = document.querySelectorAll('.code-header + .highlighter-rouge');
// const copyCodeButtons = document.querySelectorAll('.copy-code-button');

// copyCodeButtons.forEach((copyCodeButton, index) => {
//  const code = codeBlocks[index].innerText;

//  copyCodeButton.addEventListener('click', () => {
    // Copy the code to the user's clipboard
//    window.navigator.clipboard.writeText(code);

    // Update the button text visually
//    const { innerText: originalText } = copyCodeButton;
//    copyCodeButton.innerText = 'Copied!';

    // (Optional) Toggle a class for styling the button
//    copyCodeButton.classList.add('copied');

    // After 2 seconds, reset the button to its initial UI
//    setTimeout(() => {
//      copyCodeButton.innerText = originalText;
//      copyCodeButton.classList.remove('copied');
//    }, 2000);
//  });
//});

// var codeBlocks = document.querySelectorAll('pre.highlight');

// codeBlocks.forEach(function (codeBlock) {
//   var copyButton = document.createElement('button');
//   copyButton.className = 'copy';
//   copyButton.type = 'button';
//   copyButton.ariaLabel = 'Copy code to clipboard';
//   copyButton.innerText = 'Copy';

//   codeBlock.append(copyButton);

//   copyButton.addEventListener('click', function () {
//     var code = codeBlock.querySelector('code').innerText.trim();
//     window.navigator.clipboard.writeText(code);

//     copyButton.innerText = 'Copied!';
//     var fourSeconds = 2000;

//     setTimeout(function () {
//       copyButton.innerText = 'Copy';
//     }, fourSeconds);
//   });
// });


 