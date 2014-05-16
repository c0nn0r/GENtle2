/**
Handles displaying a sequence in a canvas.

Is instantiated inside a parent Backbone.View, and is automatically
rendered.

@class SequenceCanvas
@constructor
@uses SequenceCanvasContextMenu
@uses SequenceCanvasHandlers
@uses SequenceCanvasUtilities
**/ 
define(function(require) {
  'use strict';
  var Artist            = require('lib/graphics/artist'),
      Lines             = require('lib/sequence_canvas/lines'),
      Caret             = require('lib/sequence_canvas/caret'),
      CopyPasteHandler  = require('lib/copy_paste_handler'),
      _Handlers         = require('lib/sequence_canvas/_sequence_canvas_handlers'),
      _Utilities        = require('lib/sequence_canvas/_sequence_canvas_utilities'),
      _ContextMenu      = require('lib/sequence_canvas/_sequence_canvas_context_menu'),
      Hotkeys           = require('lib/hotkeys'),
      Q                 = require('q'),
      SequenceCanvas;

  SequenceCanvas = function(options) {
    var _this = this;
    options = options || {};
    this.visible = true;

    // Context binding (context is lost in Promises' `.then` and `.done`)
    _.bindAll(this, 'calculateLayoutSettings', 
                    'display', 
                    'refresh', 
                    'redraw',
                    'afterNextDisplay',
                    'handleScrolling',
                    'handleMousedown',
                    'handleMousemove',
                    'handleMouseup',
                    'handleClick',
                    'handleKeypress',
                    'handleKeydown',
                    'handleBlur'
              );

    /**
        Instance of BackboneView in which the canvas lives
        @property view
        @type Backbone.View
    **/
    this.view = options.view;

    /**
        Canvas element as a jQuery object
        @property $canvas
        @type jQuery object
        @default first CANVAS DOM element in `this.view.$el`
    **/
    this.$canvas = options.$canvas || this.view.$('canvas').first();

    /**
        Invisible DIV used to handle scrolling. 
        As high as the total virtual height of the sequence as displayed.
        @property $scrollingChild
        @type jQuery object
        @default `.scrollingChild`
    **/
    this.$scrollingChild = options.$scrollingChild || this.view.$('.scrolling-child').first();

    /**
        Div in which `this.$scrollingChild` will scroll.
        Same height as `this.$canvas`
        @property $scrollingParent
        @type jQuery object
        @default jQuery `.scrollingParent`
    **/
    this.$scrollingParent = options.$scrollingParent || this.view.$('.scrolling-parent').first();

    /**
        Sequence to be displayed
        @property sequence
        @type Sequence (Backbone.Model)
        @default `this.view.model`
    **/
    this.sequence = options.sequence || this.view.model;

    /**
        @property layoutSettings
        @type Object
        @default false
    **/
    this.layoutSettings = {
      canvasDims: {width:1138, height:448},
      pageMargins: {
        left:20,
        top:20, 
        right:20,
        bottom:20
      },
      scrollPercentage: 1.0,
      gutterWidth: 30,
      basesPerBlock: 10,
      basePairDims: {width:10, height:15},
      sequenceLength: this.sequence.length(),
      lines: options.lines || {

        // Blank line
        topSeparator: new Lines.Blank(this, {
          height: 5,
          visible: function() { return _this.sequence.get('displaySettings.rows.separators'); }
        }),

        // Position numbering
        position: new Lines.Position(this, {
          height: 15, 
          baseLine: 15, 
          textFont: "10px Monospace", 
          textColour:"#005",
          transform: function(string) {
            return string.toString().replace(/\B(?=(\d{3})+(?!\d))/g, ",");
          },
          visible: _.memoize2(function() { 
            return _this.sequence.get('displaySettings.rows.numbering'); 
          })
        }),

        // Aminoacids
        aa: new Lines.DNA(this, {
          height: 15, 
          baseLine: 15, 
          textFont: "13px Monospace", 
          transform: function(base) {
            return _this.sequence.getAA(_this.sequence.get('displaySettings.rows.aa'), base, parseInt(_this.sequence.get('displaySettings.rows.aaOffset')));
          },
          visible: _.memoize2(function() {
            return _this.sequence.get('displaySettings.rows.aa') != 'none';
          }),
          textColour: function(codon) { return {'STP': 'red', 'S  ': 'red'}[codon.sequence] || '#79B6F9'; }
        }),

        // DNA Bases
        dna: new Lines.DNA(this, {
          height: 15, 
          baseLine: 15, 
          textFont: "15px Monospace", 
          textColour:"#000",
          selectionColour: "#1a1a63",
          selectionTextColour: "#fff"
        }),

        // Complements
        complements: new Lines.DNA(this, {
          height: 15, 
          baseLine: 15, 
          textFont: "15px Monospace", 
          textColour:"#bbb",
          getSubSeq: _.partial(this.sequence.getTransformedSubSeq, 'complements', {}),
          visible: _.memoize2(function() { 
            return _this.sequence.get('displaySettings.rows.complements'); 
          })
        }),

        // Annotations
        features: new Lines.Feature(this, {
          unitHeight: 15,
          baseLine: 10,
          textFont: "11px Monospace", 
          textColour: "white",
          textPadding: 2,
          margin: 2,
          lineSize: 2,
          colour: function(type) { return {'CDS': 'blue'}[type] || 'red';},
          visible: _.memoize2(function() { 
            return _this.sequence.get('features') && _this.sequence.get('displaySettings.rows.features'); 
          })
        }),

        // Blank line
        bottomSeparator: new Lines.Blank(this, {
          height: 10,
          visible: _.memoize2(function() { 
            return _this.sequence.get('displaySettings.rows.separators'); 
          })
        })
      }
    };

    this.layoutHelpers = {};
    this.artist = new Artist(this.$canvas);
    this.caret = new Caret(this);
    this.allowedInputChars = ['A', 'T', 'C', 'G'];
    this.displayDeferred = Q.defer();
    this.copyPasteHandler = new CopyPasteHandler();

    this.contextMenu = this.view.getView('#sequence-canvas-context-menu-outlet');

    this.invertHotkeys = _.invert(Hotkeys);
    this.commandKeys = {};
    _.each(['A', 'C', 'Z', 'V'], function(key) {
      _this.commandKeys[key] = key.charCodeAt(0);
    });

    // Events
    this.view.on('resize', this.refresh);
    this.sequence.on('change:sequence', this.redraw);
    this.sequence.on('change:displaySettings.* change:features.* change:features', this.refresh);
    this.$scrollingParent.on('scroll',    this.handleScrolling);
    this.$scrollingParent.on('mousedown', this.handleMousedown);
    this.$scrollingParent.on('keypress',  this.handleKeypress);
    this.$scrollingParent.on('keydown',   this.handleKeydown);
    this.$scrollingParent.on('blur',      this.handleBlur);

    // Kickstart rendering
    this.refresh();
  };

  _.extend(SequenceCanvas.prototype, _Handlers.prototype);
  _.extend(SequenceCanvas.prototype, _Utilities.prototype);
  _.extend(SequenceCanvas.prototype, _ContextMenu.prototype);

  /**
      Updates Canvas Dimemnsions based on viewport.
      @method updateCanvasDims
      @returns {Promise} a Promise finished when this and `this.calculateLayoutSettings` are finished
  **/
  SequenceCanvas.prototype.updateCanvasDims = function() {
    var _this = this;

    return Q.promise(function(resolve, reject) {
      // Updates width of $canvas to take scrollbar of $scrollingParent into account
      _this.$canvas.width(_this.$scrollingChild.width());

      var width   = _this.$canvas[0].scrollWidth,
          height  = _this.$canvas[0].scrollHeight;

      _this.layoutSettings.canvasDims.width = width;
      _this.layoutSettings.canvasDims.height = height;

      if(_this.$canvas[0].width != width || _this.$canvas[0].height != height) {
        _this.$canvas[0].width = width;
        _this.$canvas[0].height = height;
      }

      resolve();
    });
  };


  /**
      Calculate "helper" layout settings based on already set layout settings
      @method calculateLayoutSettings
      @returns {Promise} a Promise fulfilled when finished
  **/
  SequenceCanvas.prototype.calculateLayoutSettings = function() {
    //line offsets
    var line_offset = _.values(this.layoutSettings.lines)[0].height,
        i, 
        blocks_per_row,
        ls = this.layoutSettings,
        lh = this.layoutHelpers,
        _this = this;

    return Q.promise(function(resolve, reject) {

      //basesPerRow
      blocks_per_row = Math.floor( (ls.canvasDims.width + ls.gutterWidth - (ls.pageMargins.left + ls.pageMargins.right))/(ls.basesPerBlock * ls.basePairDims.width + ls.gutterWidth) );
      if (blocks_per_row !== 0){
        lh.basesPerRow = ls.basesPerBlock*blocks_per_row;
      }else{
        lh.basesPerRow = Math.floor((ls.canvasDims.width + ls.gutterWidth - (ls.pageMargins.left + ls.pageMargins.right))/ls.basePairDims.width);
        //we want bases per row to be a multiple of 10 (DOESNT WORK)
        if (lh.basesPerRow > 5){
          lh.basesPerRow = 5;
        }else if (lh.basesPerRow > 2){
          lh.basesPerRow = 2;
        }else{
          lh.basesPerRow = 1;
        }
      }

      lh.lineOffsets = {};
      _.each(ls.lines, function(line, lineName) {
        line.clearCache();
        if(line.visible === undefined || line.visible()) {
          lh.lineOffsets[lineName] = line_offset;
          if(_.isFunction(line.calculateHeight)) line.calculateHeight();
          line_offset += line.height;
        }
      });

      //row height
      lh.rows = {height:line_offset};

      //total number of rows in sequence, 
      lh.rows.total = Math.ceil(ls.sequenceLength / lh.basesPerRow) ;
      // number of visible rows in canvas
      lh.rows.visible = Math.ceil(ls.canvasDims.height / lh.rows.height) ;

      //page dims
      lh.pageDims = {
        width:ls.canvasDims.width, 
        height: ls.pageMargins.top + ls.pageMargins.bottom + lh.rows.total*lh.rows.height 
      };

      // canvas y scrolling offset
      lh.yOffset = lh.yOffset || _this.sequence.get('displaySettings.yOffset') || 0;
      _this.$scrollingParent.scrollTop(lh.yOffset);

      _this.clearCache(); 

      // We resize `this.$scrollingChild` and fullfills the Promise
      _this.resizeScrollHelpers().then(resolve);
    });
  };

  /** 
      If `this.visible`, displays the sequence in the initiated canvas.
      @method display
  **/
  SequenceCanvas.prototype.display = function() {
    if(this.visible) {
      var context         = this.artist.context,
          ls              = this.layoutSettings,
          lh              = this.layoutHelpers,
          _this           = this,
          i, k, pos, baseRange, y;

      return Q.promise(function(resolve, reject){
        //clear canvas
        context.clearRect(0,0,context.canvas.width, context.canvas.height);

        _this.forEachRowInRange(0, ls.canvasDims.height, function(y) {
          baseRange = _this.getBaseRangeFromYPos(y);
          _.each(ls.lines, function(line, key) {
            if(line.visible === undefined || line.visible()) {
              line.draw(y, baseRange);
              y += line.height;
            }
          });
        });

        _this.displayDeferred.resolve();
        _this.displayDeferred = Q.defer();
        resolve();

      });
    } else{
      return Q.promise(function(resolve, reject) {
        reject();
      });
    }
  };

  /**
  Resizes $scrollingChild after window/parent div has been resized
  @method resizeScrollHelpers
  @return {Promise}
  **/
  SequenceCanvas.prototype.resizeScrollHelpers = function() {
    var _this = this,
        layoutHelpers = _this.layoutHelpers;
    return Q.promise(function(resolve, reject) {
      _this.$scrollingChild.height(layoutHelpers.pageDims.height);
      _this.scrollTo();
      resolve();
    });
  };

  /**
  Updates layout settings and redraws canvas
  @method refresh
  **/
  SequenceCanvas.prototype.refresh = function() {
    if(this.caretPosition) {
      this.hideCaret();
      this.caretPosition = undefined;
    }
    this.updateCanvasDims()
      .then(this.calculateLayoutSettings)
      .then(this.redraw);
  };

  /**
  Redraws canvas on the next animation frame
  @method redraw
  **/
  SequenceCanvas.prototype.redraw = function() {
    return requestAnimationFrame(this.display);
  };

  SequenceCanvas.prototype.scrollTo = function(yOffset) {
    var deferred = Q.defer();

    if(yOffset !== undefined) {
      this.sequence.set('displaySettings.yOffset', 
        this.layoutHelpers.yOffset = yOffset,
        { silent: true }
      );
      this.sequence.throttledSave();
    }

    this.$scrollingParent.scrollTop(this.layoutHelpers.yOffset);

    this.afterNextDisplay(deferred.resolve);

    this.redraw();

    this.restoreContextMenuYPosition();

    return deferred.promise;
  };

  /**
  Make base visible (if it is below the visible part of the canvas,
  will just scroll down one row)
  @method scrollBaseToVisibility
  **/
  SequenceCanvas.prototype.scrollBaseToVisibility = function(base) {
    var distanceToVisibleCanvas = this.distanceToVisibleCanvas(base);

    if(distanceToVisibleCanvas !== 0) {
      return this.scrollTo(this.layoutHelpers.yOffset + distanceToVisibleCanvas);
    } else {
      return Q.resolve();
    }
  };

  SequenceCanvas.prototype.scrollToBase = function(base) {
    if(!this.isBaseVisible(base)) {
      var yPos = this.getYPosFromBase(base),
          maxY = this.$scrollingChild.height() - this.$scrollingParent.height();
      return this.scrollTo(Math.min(yPos, maxY));
    } else {
      return Q.resolve();
    }
  };

  SequenceCanvas.prototype.clearCache = function() {
    this.getXPosFromBase.cache = {};
    // this.getYPosFromBase.cache = {};
  };

  SequenceCanvas.prototype.afterNextDisplay = function() {
    var _this = this,
        args = _.toArray(arguments),
        func = args.shift();

    this.displayDeferred.promise.then(function() {
      func.apply(_this, args);
    });
  };

  /**
  Displays the caret before a base
  @method displayCaret
  @param base [base] 
  **/
  SequenceCanvas.prototype.displayCaret = function(base) {
    var layoutHelpers = this.layoutHelpers,
        lineOffsets   = layoutHelpers.lineOffsets,
        yOffset       = layoutHelpers.yOffset,
        _this         = this,
        posX, posY;

    if(base === undefined && this.caretPosition) {
      base = this.caretPosition;
    }

    if(base > this.sequence.length()) {
      base = this.sequence.length();
    }

    this.scrollBaseToVisibility(base).then(function() {

      posX = _this.getXPosFromBase(base);
      posY = _this.getYPosFromBase(base) + lineOffsets.dna;

      _this.caret.move(posX, posY);
      _this.caretPosition = base;
      _this.showContextMenuButton(posX, posY + 20);

    });
  
  };

  SequenceCanvas.prototype.displayCaretAfterNextDisplay = 
    _.wrap(
      SequenceCanvas.prototype.displayCaret,
      SequenceCanvas.prototype.afterNextDisplay
    );

  SequenceCanvas.prototype.hideCaret = function(hideContextMenu) {
    this.caret.remove();
    if(hideContextMenu === true) {
      this.hideContextMenuButton();
    }
  };

  /**
  @method select
  **/
  SequenceCanvas.prototype.select = function(start, end) {
    this.hideCaret();
    if(start !== undefined) {
      if(start < end) {
        this.selection = [start, end];
        this.caretPosition = end + 1;
      } else {
        this.selection = [end, start];
        this.caretPosition = start + 1;
      }
    } else {
      this.selection = undefined;
      this.caretPosition = undefined;
    }
    this.redraw();
  };

  SequenceCanvas.prototype.expandSelectionToNewCaret = function(newCaret) {
    var selection = this.selection,
        previousCaret = this.caretPosition;

    if(selection[0] == selection[1] && (
        (previousCaret > selection[0] && newCaret == selection[0]) ||
        (previousCaret == selection[0] && newCaret == selection[0] + 1)
      )) {
      this.select(undefined);
    } else {
      if(newCaret > selection[0]) {
        if(previousCaret <= selection[0]) {
          this.select(newCaret, selection[1]);
        } else {
          this.select(selection[0], newCaret - 1);
        }
      } else {
        if(previousCaret <= selection[1] && newCaret < selection[1]) {
          this.select(newCaret, selection[1]);
        } else {
          this.select(newCaret, selection[0] - 1);
        }
      }
    }
    this.displayCaretAfterNextDisplay(newCaret);
  };

  SequenceCanvas.prototype.cleanPastedText = function(text) {
    var regexp = new RegExp('[^' + this.allowedInputChars.join('') + ']', 'g')
    return text.toUpperCase().replace(regexp, '');
  };

  SequenceCanvas.prototype.focus = function() {
    this.$scrollingParent.focus();
  };




  

  return SequenceCanvas;
});