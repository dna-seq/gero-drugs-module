widgetGenerators['pharmgkb'] = {
	'variant': {
		'width': 680, 
		'height': 280, 
		'function': function (div, row, tabName) {
			var chemical = getWidgetData(tabName, 'pharmgkb', row, 'chemical');
			var assocs = getWidgetData(tabName, 'pharmgkb', row, 'drug_assoc');
			if (assocs != undefined && assocs != null) {
				var pharmId = getWidgetData(tabName, 'pharmgkb', row, 'id');
				addInfoLineLink(div, 'Variant', pharmId, `https://pharmgkb.org/variant/${pharmId}`);
				assocs = JSON.parse(assocs);
				var table = getWidgetTableFrame();
				addEl(div, table);
				table.style.tableLayout = 'auto';
				table.style.width = '100%';
				var thead = getWidgetTableHead(
					['Chemicals', 'Description','Category','Significant','Study', 'Notes'], 
					['15%',       '30%',     '8%',     '10%',        '12%',   '25%'  ]
				);
				addEl(table, thead);
				var tbody = getEl('tbody');
				addEl(table, tbody);
				for (let row of assocs) {
					let pmLink = getEl('a');
					pmLink.href = `https://pubmed.ncbi.nlm.nih.gov/${row[4]}`;
					pmLink.innerText = row[4];
					pmLink.target = '_blank';
					let chemSpan = getEl('span');
					for (let i=0; i<row[0].length; i++) {
						let chemInfo = row[0][i];
						let chemLink = getEl('a');
						chemLink.innerText = chemInfo[0];
						chemLink.href = chemInfo[1];
						chemLink.target = '_blank';
						addEl(chemSpan, chemLink);
						if (i === row[0].length-2) {
							addEl(chemSpan, getTn(', '));
						}
					}
					var tr = elementsToWidgetTr([chemSpan, row[1], row[2], row[3], pmLink, row[5]]);
					addEl(tbody, tr);
				}
				addEl(div, addEl(table, tbody));
			} else if (chemical != undefined && chemical != null && chemical.indexOf('[[') == 0) {
				var chemsDiv = getEl('div');
				addEl(div, chemsDiv);
				var titleSpan = getEl('span');
				addEl(chemsDiv, titleSpan);
				addEl(titleSpan, getTn('Chemicals: '));
				titleSpan.style['font-weight'] = 'bold';
				var chemicals = JSON.parse(chemical);
				for (let i=0; i<chemicals.length; i++) {
                    var row = chemicals[i];
					var curchem = row[0];
					var curid = row[1];
					var url = `https://pharmgkb.org/chemical/${curid}`;
					var link = getEl('a');
					addEl(chemsDiv, link);
					link.href = url;
					link.target = '_blank';
					addEl(link, getTn(curchem));
					if (i + 1 < chemicals.length) {
						addEl(chemsDiv, getTn(', '));
					}
				}
				addInfoLine(div, 'Description', getWidgetData(tabName, 'pharmgkb', row, 'sentence'), tabName);
			} else {
                var chemical = getWidgetData(tabName, 'pharmgkb', row, 'chemical');
                var chemid = getWidgetData(tabName, 'pharmgkb', row, 'chemid');
                var chemsDiv = getEl('div');
                addEl(div, chemsDiv);
                var titleSpan = getEl('span');
                addEl(chemsDiv, titleSpan);
                addEl(titleSpan, getTn('Chemicals: '));
                titleSpan.style['font-weight'] = 'bold';
                if (chemical && chemid) {
                    var chemicals = chemical.split(';');
                    var chemids = chemid.split(';');
                    for (let i=0; i<chemicals.length; i++) {
                        var curchem = chemicals[i];
                        var curid = chemids[i];
                        var url = `https://pharmgkb.org/chemical/${curid}`;
                        var link = getEl('a');
                        addEl(chemsDiv, link);
                        link.href = url;
                        link.target = '_blank';
                        addEl(link, getTn(curchem));
                        if (i + 1 < chemicals.length) {
                            addEl(chemsDiv, getTn(', '));
                        }
                    }
                }
				addInfoLine(div, 'Description', getWidgetData(tabName, 'pharmgkb', row, 'sentence'), tabName);
            }
		}
	}
}
