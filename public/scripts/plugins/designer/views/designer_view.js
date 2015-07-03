import Backbone from 'backbone';
import _ from 'underscore';
import AssembleSequence from '../lib/assemble_sequence';
import template from '../templates/designer_view_template.hbs';
import AvailableSequencesView from './available_sequences_view';
import DesignedSequenceView from './designed_sequence_view';
import Gentle from 'gentle';
import uploadMultipleSequences from '../lib/upload_multiple_sequences';
import Modal from '../../../common/views/modal_view';
import DiagnosticModalView from './designer_diagnostic_modal_view';

var filterSequencesByStickyEnds = function(assembleSequence, [stickyEndStartName, stickyEndEndName]) {
  var stickyEndName = `${stickyEndStartName}-${stickyEndEndName}`.toLowerCase();
  return function() {
    return _.filter(assembleSequence.allSequences, function(sequence) {
      var stickyEnds = sequence.getStickyEnds();
      return stickyEnds && 
        `${stickyEnds.start.name}-${stickyEnds.end.name}`.toLowerCase() === stickyEndName;
    });
  };
};

var DesignerView = Backbone.View.extend({
  template: template,
  manage: true,
  className: 'designer',

  events: {
    // 'click .toggle-annotations': 'toggleAnnotations',
    // 'click .toggle-uninsertable-sequences': 'toggleUninsertableSequences',
    'change #circularise-dna': 'updateCirculariseDna',
    'click .designer-available-sequences-header button': 'triggerFileInput',
    'change .file-upload-input': 'uploadNewSequences',
    'click .assemble-sequence-btn': 'assembleSequence'
  },

  initialize: function() {
    this.model = new AssembleSequence(Gentle.currentSequence);

    // Default to `circular` true
    if(this.model.sequences.length === 0) {
      this.model.set({'isCircular': true}, {silent: true});
    }

    this.setView(
      '.designer-available-sequences-outlet.outlet-1', 
      new AvailableSequencesView({
        name: 'x-z\' Parts',
        getSequences: filterSequencesByStickyEnds(this.model, ['x', 'z\''])
      })
    );

    this.setView(
      '.designer-available-sequences-outlet.outlet-2', 
      new AvailableSequencesView({
        name: 'z-x\' Parts',
        getSequences: filterSequencesByStickyEnds(this.model, ['z', 'x\''])
      })
    );

    var designedSequenceView = this.designedSequenceView = 
      new DesignedSequenceView({model: this.model});
    this.setView('.designer-designed-sequence-outlet', designedSequenceView);
  },

  triggerFileInput: function(event) {
    event.preventDefault();
    this.$('.file-upload-input').click();
  },

  uploadNewSequences: function(event) {
    uploadMultipleSequences(event.target.files).then((sequences) => {
      this.model.addSequences(sequences);
      this.render();
    }).done();
  },

  serialize: function() {
    return {
      sequenceName: this.model.get('name'),
      submitDisabled: this.model.sequences.length === 0,
      circulariseDna: this.model.get('isCircular'),
      // insertableSequences: _.pluck(this.model.insertableSequences, 'id'),
      // uninsertableSequences: this.model.incompatibleSequences.length + this.model.lackStickyEndSequences.length,
      // incompatibleSequences: _.pluck(this.model.incompatibleSequences, 'id'),
      // lackStickyEndSequences: _.pluck(this.model.lackStickyEndSequences, 'id'),
      // showAnnotations: Gentle.currentUser.get('displaySettings.designerView.showAnnotations') || false,
      // showUninsertableSequences: Gentle.currentUser.get('displaySettings.designerView.showUninsertableSequences') || false,
    };
  },

  updateCirculariseDna: function(event) {
    event.preventDefault();
    this.model.set('isCircular', event.target.checked).throttledSave();
  },

  insertSequenceViews: function() {
    var _this = this,
        designedSequenceView;

    _.each(this.model.allSequences, function(sequence) {
      var outletSelector = `.designer-available-sequence-outlet[data-sequence_id="${sequence.id}"]`;
      var sequenceView = new AvailableSequenceView({model: sequence});
      _this.setView(outletSelector, sequenceView);
      sequenceView.render();
    });

    designedSequenceView = new DesignedSequenceView({model: this.model});
    this.setView('.designer-designed-sequence-outlet', designedSequenceView);
    this.designedSequenceView = designedSequenceView;
    designedSequenceView.render();
  }, 

  getAvailableSequenceViewFromSequenceId: function(sequenceId) {
    return this.getView(`.designer-available-sequence-outlet[data-sequence_id="${sequenceId}"]`);
  },

  hoveredOverSequence: function(sequenceId) {
    var indices = this.model.insertabilityState[sequenceId];
    this.designedSequenceView.highlightDropSites(indices);
  },

  unhoveredOverSequence: function(sequenceId) {
    this.designedSequenceView.unhighlightDropSites();
  },

  // toggleAnnotations: function(event) {
  //   var showAnnotations = Gentle.currentUser.get('displaySettings.designerView.showAnnotations');
  //   showAnnotations = _.isUndefined(showAnnotations) ? true : !showAnnotations;
  //   Gentle.currentUser.set('displaySettings.designerView.showAnnotations', showAnnotations);
  //   this.render();
  // },

  // toggleUninsertableSequences: function(event) {
  //   var showUninsertableSequences = Gentle.currentUser.get('displaySettings.designerView.showUninsertableSequences');
  //   showUninsertableSequences = _.isUndefined(showUninsertableSequences) ? true : !showUninsertableSequences;
  //   Gentle.currentUser.set('displaySettings.designerView.showUninsertableSequences', showUninsertableSequences);
  //   this.render();
  // },

  isInsertable: function(sequence) {
    return this.model.isInsertable(sequence);
  },

  getDescriptiveAnnotationContent: function(sequence) {
    var features = sequence.get('features');
    if(features.length == 1) {
      var feature = features[0];
      var range = feature.ranges[0];
      if(range.from === 0 && range.to >= sequence.length()-1) {
        return feature.name;
      }
    }
  },

  changeSecondaryView: function() {
    // Currently NoOp
  },

  cleanup: function() {
    this.removeAllViews();
    Gentle.sequences.off(null, null, this);
  },

  removeAllViews: function() {
    this.designedSequenceView = undefined;
    this.getViews().each((view) => {
      view.remove();
    });
  },

  assembleSequence: function() {
    var errors = this.model.errors;
    if(errors.length > 0) {
      Modal.show({
        title: 'Cannot assemble circuit',
        confirmLabel: 'OK',
        cancelLabel: null,
        bodyView: new DiagnosticModalView({
          errors: _.uniq(_.pluck(errors, 'type'))
        })
      });
    }
  }

});

export default DesignerView;